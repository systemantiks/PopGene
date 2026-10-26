/* PopGene's bounded ABI around the source-vendored Primer3 library.
 * All calls are serialized by primer.rs because Primer3 uses global caches. */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "libprimer3.h"
#include "oligotm.h"
#include "thal.h"
#include "thal_parameters.h"

typedef struct {
    double dna_nm, monovalent_mm, magnesium_mm, dntp_mm, dmso_percent,
           formamide_m, temperature_c;
    int tm_method, salt_method;
} pg_reaction;
typedef struct {
    double tm, dh, ds, dg;
    int found;
} pg_structure;
typedef struct {
    pg_reaction reaction;
    double min_tm, optimal_tm, max_tm, min_gc, max_gc, max_tm_difference;
    int min_length, optimal_length, max_length, product_min, product_max,
        num_return, target_start, target_length;
} pg_design_settings;
typedef struct {
    int forward_start, forward_length, reverse_end, reverse_length, product_size;
    double forward_tm, reverse_tm, forward_gc, reverse_gc, penalty,
           forward_hairpin, reverse_hairpin, pair_any, pair_end;
} pg_pair;

static int pg_initialize(char *error, size_t capacity) {
    static int initialized = 0;
    if (initialized) return 0;
    thal_parameters parameters;
    thal_results result = {0};
    thal_set_null_parameters(&parameters);
    if (set_default_thal_parameters(&parameters) != 0) {
        snprintf(error, capacity, "Cannot initialize Primer3 thermodynamic parameters");
        thal_free_parameters(&parameters);
        return -1;
    }
    int status = get_thermodynamic_values(&parameters, &result);
    thal_free_parameters(&parameters);
    if (status != 0) { snprintf(error, capacity, "%s", result.msg); return -1; }
    initialized = 1;
    return 0;
}

double pg_primer_tm(const char *sequence, const pg_reaction *r) {
    return oligotm(sequence, r->dna_nm, r->monovalent_mm, r->magnesium_mm,
                   r->dntp_mm, r->dmso_percent, 0.6, r->formamide_m,
                   (tm_method_type)r->tm_method,
                   (salt_correction_type)r->salt_method, -10.0).Tm;
}

int pg_primer_structure(const char *first, const char *second,
                        const pg_reaction *r, int hairpin,
                        pg_structure *output, char *error, size_t capacity) {
    if (pg_initialize(error, capacity)) return -1;
    thal_args args;
    thal_results result = {0};
    set_thal_default_args(&args);
    args.type = hairpin ? thal_hairpin : thal_any;
    args.mv = r->monovalent_mm;
    args.dv = r->magnesium_mm;
    args.dntp = r->dntp_mm;
    args.dna_conc = r->dna_nm;
    args.temp = r->temperature_c + 273.15;
    thal((const unsigned char *)first, (const unsigned char *)second,
         &args, THL_GENERAL, &result, 0);
    if (result.temp == THAL_ERROR_SCORE || !isfinite(result.temp)) {
        snprintf(error, capacity, "%s", result.msg); return -1;
    }
    output->tm = result.temp;
    output->dh = result.dh;
    output->ds = result.ds;
    output->dg = result.dg;
    output->found = !result.no_structure;
    return 0;
}

int pg_design_primers(const char *sequence, const pg_design_settings *s,
                     pg_pair *pairs, int capacity, int *count,
                     char *error, size_t error_capacity,
                     char *warning, size_t warning_capacity) {
    *count = 0;
    if (pg_initialize(error, error_capacity)) return -1;
    p3_global_settings *settings = p3_create_global_settings();
    seq_args_t *args = create_seq_arg();
    p3retval *result = NULL;
    int status = -1;
    if (!settings || !args) {
        snprintf(error, error_capacity, "Primer3 allocation failed"); goto cleanup;
    }
    settings->p_args.min_size = s->min_length;
    settings->p_args.opt_size = s->optimal_length;
    settings->p_args.max_size = s->max_length;
    settings->p_args.min_tm = s->min_tm;
    settings->p_args.opt_tm = s->optimal_tm;
    settings->p_args.max_tm = s->max_tm;
    settings->p_args.min_gc = s->min_gc;
    settings->p_args.max_gc = s->max_gc;
    settings->p_args.salt_conc = s->reaction.monovalent_mm;
    settings->p_args.divalent_conc = s->reaction.magnesium_mm;
    settings->p_args.dntp_conc = s->reaction.dntp_mm;
    settings->p_args.dna_conc = s->reaction.dna_nm;
    settings->p_args.dmso_conc = s->reaction.dmso_percent;
    settings->p_args.dmso_fact = 0.6;
    settings->p_args.formamide_conc = s->reaction.formamide_m;
    settings->tm_santalucia = (tm_method_type)s->reaction.tm_method;
    settings->salt_corrections = (salt_correction_type)s->reaction.salt_method;
    settings->max_diff_tm = s->max_tm_difference;
    settings->num_return = s->num_return;
    settings->thermodynamic_oligo_alignment = 1;
    settings->thermodynamic_template_alignment = 0;
    p3_empty_gs_product_size_range(settings);
    p3_add_to_gs_product_size_range(settings, s->product_min, s->product_max);
    if (p3_set_sa_sequence(args, sequence)) {
        snprintf(error, error_capacity, "Primer3 sequence allocation failed"); goto cleanup;
    }
    if (s->target_length && p3_add_to_sa_tar2(args, s->target_start, s->target_length)) {
        snprintf(error, error_capacity, "Invalid target interval"); goto cleanup;
    }
    result = choose_primers(settings, args);
    if (!result) { snprintf(error, error_capacity, "Primer3 design failed"); goto cleanup; }
    const char *global = p3_get_rv_global_errors(result);
    const char *local = p3_get_rv_per_sequence_errors(result);
    if ((global && *global) || (local && *local)) {
        snprintf(error, error_capacity, "%s%s%s", global ? global : "",
                 global && local ? "; " : "", local ? local : ""); goto cleanup;
    }
    const char *warnings = p3_get_rv_warnings(result);
    if (warnings) snprintf(warning, warning_capacity, "%s", warnings);
    const pair_array_t *best = p3_get_rv_best_pairs(result);
    for (int i = 0; i < best->num_pairs && i < capacity; ++i) {
        const primer_pair *p = &best->pairs[i];
        pairs[i] = (pg_pair) {
            p->left->start, p->left->length, p->right->start, p->right->length,
            p->product_size, p->left->temp, p->right->temp,
            p->left->gc_content, p->right->gc_content, p->pair_quality,
            p->left->hairpin_th, p->right->hairpin_th, p->compl_any, p->compl_end
        };
        ++*count;
    }
    status = 0;
cleanup:
    if (result) { destroy_secundary_structures(settings, result); destroy_p3retval(result); }
    destroy_dpal_thal_arg_holder();
    if (args) destroy_seq_args(args);
    if (settings) p3_destroy_global_settings(settings);
    return status;
}
