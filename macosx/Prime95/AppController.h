//
//  AppController.h
//  Prime95
//
//  Created by George Woltman on 4/17/09.
//  Copyright 2009-2017 Mersenne Research, Inc. All rights reserved.
//

#import <Cocoa/Cocoa.h>
@class AboutController;
@class BenchmarkController;
@class ContinueController;
@class CPUController;
@class ECMController;
@class ManualCommunicationController;
@class Pminus1Controller;
@class PreferencesController;
@class PrimeNetController;
@class PRPController;
@class QuitGIMPSController;
@class StatusController;
@class StopController;
@class TestController;
@class TimeController;
@class TortureTestController;
@class UnreserveController;
@class WorkerWindowsController;

@interface AppController : NSObject {
	AboutController *aboutController;
	BenchmarkController *benchmarkController;
	ContinueController *continueController;
	CPUController *cpuController;
	ECMController *ecmController;
	ManualCommunicationController *manualCommunicationController;
	Pminus1Controller *pminus1Controller;
	PreferencesController *preferencesController;
	PrimeNetController *primeNetController;
	PRPController *prpController;
	QuitGIMPSController *quitGIMPSController;
	StatusController *statusController;
	StopController *stopController;
	TestController *testController;
	TimeController *timeController;
	TortureTestController *tortureTestController;
	UnreserveController *unreserveController;
	WorkerWindowsController *workerWindowsController;
}

- (IBAction)showAboutPanel:(id)sender;
- (IBAction)showAboutPrimenet:(id)sender;
- (IBAction)shutDown:(id)sender;
- (IBAction)testPrimeNet:(id)sender;
- (IBAction)testWorkerWindows:(id)sender;
- (IBAction)testStatus:(id)sender;
- (IBAction)testContinue:(id)sender;
- (IBAction)testStop:(id)sender;
- (IBAction)editBigger:(id)sender;
- (IBAction)editSmaller:(id)sender;
- (IBAction)advancedTest:(id)sender;
- (IBAction)advancedTime:(id)sender;
- (IBAction)advancedPminus1:(id)sender;
- (IBAction)advancedECM:(id)sender;
- (IBAction)advancedPRP:(id)sender;
- (IBAction)advancedManualCommunication:(id)sender;
- (IBAction)toggleSuminpErrorChecking:(id)sender;
- (IBAction)toggleErrorChecking:(id)sender;
- (IBAction)advancedUnreserve:(id)sender;
- (IBAction)advancedQuitGIMPS:(id)sender;
- (IBAction)optionsCPU:(id)sender;
- (IBAction)optionsPreferences:(id)sender;
- (IBAction)optionsTortureTest:(id)sender;
- (IBAction)optionsBenchmark:(id)sender;
- (IBAction)toggleMergeMainComm:(id)sender;
- (IBAction)toggleMergeMainCommWorker:(id)sender;
- (IBAction)toggleMergeAllWorkers:(id)sender;
- (IBAction)helpMersenneForum:(id)sender;
- (IBAction)helpMersenneWiki:(id)sender;

@end

extern AppController *myAppController;

