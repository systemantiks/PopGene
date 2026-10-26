// Sparkle is loaded only from our signed app bundle. Keeping the bridge dynamically
// linked lets `cargo run/test` work without downloading or installing Sparkle.
#import <AppKit/AppKit.h>
#import <objc/runtime.h>

// The small runtime interface we use, matching Sparkle 2.10's public headers.
// No Sparkle class symbols are referenced at link time.
@protocol PGUpdater <NSObject>
- (BOOL)startUpdater:(NSError **)error;
@property BOOL automaticallyChecksForUpdates;
@property BOOL automaticallyDownloadsUpdates;
@end
@protocol PGUpdaterController <NSObject>
- (instancetype)initWithStartingUpdater:(BOOL)start
                       updaterDelegate:(id)delegate
                    userDriverDelegate:(id)driverDelegate;
@property(readonly) id<PGUpdater> updater;
- (void)checkForUpdates:(id)sender;
@end

static id<PGUpdaterController> controller;
static NSString *lastError;
static BOOL quitRequested;
static BOOL quitCancelled;
static BOOL quitApproved;
static BOOL quitGuardInstalled;
static void (^pendingInstall)(void);

@interface PGUpdateDelegate : NSObject
@end
@implementation PGUpdateDelegate
- (void)updater:(id)updater didAbortWithError:(NSError *)error {
    (void)updater;
    (void)error;
    if (!pendingInstall && !quitApproved) return;
    pendingInstall = nil;
    quitApproved = NO;
    quitRequested = NO;
    quitCancelled = YES;
}
- (BOOL)updater:(id)updater shouldPostponeRelaunchForUpdate:(id)item
        untilInvokingBlock:(void (^)(void))installHandler {
    (void)updater;
    (void)item;
    pendingInstall = [installHandler copy];
    quitRequested = YES;
    return YES;
}
@end
static PGUpdateDelegate *updateDelegate;

// Add the standard AppKit termination hook to GPUI's delegate without replacing
// it: GPUI retains all its menu, window, lifecycle and file-open callbacks.
// Sparkle terminates via NSApp, so this covers every updater installation path,
// including retries after cancellation and an update already staged on quit.
static NSApplicationTerminateReply shouldTerminate(id self, SEL command, NSApplication *app) {
    (void)self;
    (void)command;
    (void)app;
    if (quitApproved) return NSTerminateNow;
    quitRequested = YES;
    // NSTerminateLater enters AppKit's modal quit loop, which prevents GPUI's
    // dispatch executor from rendering its save dialog. Cancel this attempt and
    // let GPUI retry once it has saved/discarded. Sparkle normally reaches us only
    // after the postponement delegate above has completed that process.
    return NSTerminateCancel;
}

int popgene_install_quit_guard(void) {
    NSCAssert(NSThread.isMainThread, @"Updater must run on the main thread");
    if (quitGuardInstalled) return 1;
    id delegate = NSApp.delegate;
    if (!delegate || [delegate respondsToSelector:@selector(applicationShouldTerminate:)]) return 0;
    // Match the platform ABI rather than hard-coding NSInteger's encoding.
    NSString *types = [NSString stringWithFormat:@"%s@:@", @encode(NSApplicationTerminateReply)];
    quitGuardInstalled = class_addMethod([delegate class], @selector(applicationShouldTerminate:),
                                        (IMP)shouldTerminate, types.UTF8String);
    return quitGuardInstalled;
}

int popgene_take_quit_request(void) {
    if (quitCancelled) {
        quitCancelled = NO;
        return 2;
    }
    BOOL requested = quitRequested;
    quitRequested = NO;
    return requested;
}

// A cancelled restart retains Sparkle's installation block. A later manual
// update check or normal Quit retries it. Never call it until all documents are
// saved or explicitly discarded. Dispatch outside the current GPUI borrow.
int popgene_reply_to_quit(int approved) {
    quitApproved = approved != 0;
    quitRequested = NO;
    if (!approved || !pendingInstall) return 0;
    void (^install)(void) = pendingInstall;
    pendingInstall = nil;
    dispatch_async(dispatch_get_main_queue(), install);
    return 1;
}

const char *popgene_check_for_updates(void) {
    NSCAssert(NSThread.isMainThread, @"Updater must run on the main thread");
    if (!quitGuardInstalled) {
        return "Updates are unavailable because PopGene could not protect unsaved documents during restart.";
    }
    if (pendingInstall) {
        quitRequested = YES;
        return NULL;
    }
    if (!controller) {
        NSBundle *app = NSBundle.mainBundle;
        NSString *feed = [app objectForInfoDictionaryKey:@"SUFeedURL"];
        NSString *key = [app objectForInfoDictionaryKey:@"SUPublicEDKey"];
        if (!feed.length || !key.length) {
            return "This development build is not configured for updates. Install a published PopGene release to receive updates.";
        }
        NSURL *frameworkURL = [app.privateFrameworksURL URLByAppendingPathComponent:@"Sparkle.framework"];
        NSBundle *framework = [NSBundle bundleWithURL:frameworkURL];
        NSError *error = nil;
        if (!framework || ![framework loadAndReturnError:&error]) {
            lastError = [NSString stringWithFormat:@"Could not load the update framework: %@",
                         error.localizedDescription ?: @"Sparkle.framework is missing from this app."];
            return lastError.UTF8String;
        }
        Class updaterClass = NSClassFromString(@"SPUStandardUpdaterController");
        if (!updaterClass) return "The bundled update framework is incompatible with PopGene.";
        updateDelegate = [PGUpdateDelegate new];
        controller = [(id<PGUpdaterController>)[updaterClass alloc]
                      initWithStartingUpdater:NO updaterDelegate:updateDelegate userDriverDelegate:nil];
        // Manual checks only, including when old preferences exist on this Mac.
        controller.updater.automaticallyChecksForUpdates = NO;
        controller.updater.automaticallyDownloadsUpdates = NO;
        if (![controller.updater startUpdater:&error]) {
            controller = nil;
            lastError = [NSString stringWithFormat:@"Could not start the updater: %@", error.localizedDescription];
            return lastError.UTF8String;
        }
    }
    [controller checkForUpdates:nil];
    return NULL;
}
