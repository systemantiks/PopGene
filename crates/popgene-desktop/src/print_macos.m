// AppKit print sheet for PDF pages rendered by the figure pipeline.
#import <AppKit/AppKit.h>
#import <PDFKit/PDFKit.h>

static BOOL printing;
// The sheet runs after popgene_print_pdf returns; keep the pages alive until it finishes.
static PDFDocument *printDocument;

@interface PGPrintDelegate : NSObject
@end
@implementation PGPrintDelegate
- (void)printOperationDidRun:(NSPrintOperation *)operation
                     success:(BOOL)success
                 contextInfo:(void *)info {
    (void)operation;
    (void)success;
    (void)info;
    printDocument = nil;
    printing = NO;
}
@end
static PGPrintDelegate *printDelegate;

// Copies the PDF bytes. Presentation is queued on the main queue so AppKit never
// runs the sheet while GPUI is updating its entities.
const char *popgene_print_pdf(const unsigned char *bytes, size_t length, const char *title) {
    NSCAssert(NSThread.isMainThread, @"Printing must start on the main thread");
    if (printing) return "A print dialog is already open.";
    PDFDocument *document =
        [[PDFDocument alloc] initWithData:[NSData dataWithBytes:bytes length:length]];
    if (!document || document.pageCount == 0) return "The print pages could not be prepared.";
    NSString *jobTitle = title ? @(title) : nil;
    printing = YES;
    dispatch_async(dispatch_get_main_queue(), ^{
        NSPrintInfo *info = [NSPrintInfo.sharedPrintInfo copy];
        NSPrintOperation *operation = [document printOperationForPrintInfo:info
                                                               scalingMode:kPDFPrintPageScaleDownToFit
                                                                autoRotate:YES];
        if (!operation) {
            printing = NO;
            return;
        }
        if (jobTitle.length) operation.jobTitle = jobTitle;
        operation.showsPrintPanel = YES;
        operation.showsProgressPanel = YES;
        operation.printPanel.options |= NSPrintPanelShowsPaperSize | NSPrintPanelShowsOrientation |
                                        NSPrintPanelShowsScaling | NSPrintPanelShowsPreview;
        printDocument = document;
        if (!printDelegate) printDelegate = [PGPrintDelegate new];
        NSWindow *window =
            NSApp.keyWindow ?: NSApp.mainWindow ?: NSApp.orderedWindows.firstObject;
        if (window) {
            [operation runOperationModalForWindow:window
                                         delegate:printDelegate
                                   didRunSelector:@selector(printOperationDidRun:success:contextInfo:)
                                      contextInfo:NULL];
        } else {
            // No window can host a sheet; this block is outside GPUI's update cycle.
            [operation runOperation];
            printDocument = nil;
            printing = NO;
        }
    });
    return NULL;
}
