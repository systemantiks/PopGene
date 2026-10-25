// Integration test of the actual AppKit termination protocol used by Sparkle.
// Build with updater_macos.m, run in a graphical macOS login session. No updates
// are downloaded and no application is replaced by this test.
#import <AppKit/AppKit.h>
#include <assert.h>

extern int popgene_install_quit_guard(void);
extern int popgene_take_quit_request(void);
extern int popgene_reply_to_quit(int approved);
extern const char *popgene_check_for_updates(void);

@interface PGUpdateDelegate : NSObject
- (BOOL)updater:(id)updater shouldPostponeRelaunchForUpdate:(id)item untilInvokingBlock:(void (^)(void))block;
- (void)updater:(id)updater didAbortWithError:(NSError *)error;
@end

@interface TerminationTest : NSObject <NSApplicationDelegate>
@end
@implementation TerminationTest
- (void)applicationDidFinishLaunching:(NSNotification *)notification {
    (void)notification;
    assert(popgene_install_quit_guard());
    assert(popgene_check_for_updates() != NULL); // Unbundled development build.
    PGUpdateDelegate *failedUpdater = [PGUpdateDelegate new];
    [failedUpdater updater:nil didAbortWithError:nil];
    assert(popgene_take_quit_request() == 0); // A check error does not cancel editing.
    assert([failedUpdater updater:nil shouldPostponeRelaunchForUpdate:nil untilInvokingBlock:^{
        abort(); // An aborted install must never resume this block.
    }]);
    [failedUpdater updater:nil didAbortWithError:nil];
    assert(popgene_take_quit_request() == 2); // Reset a previously approved close.
    assert(popgene_take_quit_request() == 0);
    dispatch_async(dispatch_get_main_queue(), ^{ [NSApp terminate:nil]; });
    dispatch_after(dispatch_time(DISPATCH_TIME_NOW, NSEC_PER_SEC / 4), dispatch_get_main_queue(), ^{
        assert(popgene_take_quit_request());
        assert(!popgene_take_quit_request());
        assert(!popgene_reply_to_quit(0)); // Cancel must leave the app running.
        PGUpdateDelegate *updater = [PGUpdateDelegate new];
        assert([updater updater:nil shouldPostponeRelaunchForUpdate:nil untilInvokingBlock:^{
            [NSApp terminate:nil];
        }]);
        assert(popgene_take_quit_request());
        assert(!popgene_reply_to_quit(0)); // Cancel retains the pending install.
        dispatch_after(dispatch_time(DISPATCH_TIME_NOW, NSEC_PER_SEC / 4), dispatch_get_main_queue(), ^{
            assert(!popgene_take_quit_request());
            assert(popgene_check_for_updates() == NULL); // Retry the staged update.
        });
        dispatch_after(dispatch_time(DISPATCH_TIME_NOW, NSEC_PER_SEC / 2), dispatch_get_main_queue(), ^{
            assert(popgene_take_quit_request()); // Sparkle can retry a cancelled exit.
            assert(popgene_reply_to_quit(1)); // Resume its block, then terminate.
        });
    });
}
- (void)applicationWillTerminate:(NSNotification *)notification {
    (void)notification;
    puts("PASS: native termination deferred, cancelled, retried and approved");
}
@end

int main(void) {
    @autoreleasepool {
        [NSApplication sharedApplication];
        [NSApp setActivationPolicy:NSApplicationActivationPolicyAccessory];
        TerminationTest *delegate = [TerminationTest new];
        NSApp.delegate = delegate;
        dispatch_after(dispatch_time(DISPATCH_TIME_NOW, 5 * NSEC_PER_SEC), dispatch_get_main_queue(), ^{
            fputs("FAIL: termination test timed out\n", stderr);
            exit(1);
        });
        [NSApp run];
    }
    return 0;
}
