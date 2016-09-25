//
//  BundlerTrackViewerAppDelegate.h
//  BundlerTrackViewer
//
//  Created by Duygu Ceylan on 11/8/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import "OpenGLViewPort.h"

class BundlerOutputParser;
class BundlerManager;

enum VISUALIZATION_TYPE
{
	FeaturePoints = 0,
	Repetitions,
	FeatureMatches,
	Tracks,
	GridStrokes
};

@interface BundlerTrackViewerAppDelegate : NSObject <NSApplicationDelegate> {
    NSWindow *window;
	
	BundlerManager *bundlerManager;
	BundlerOutputParser *bundleParser;
	
	IBOutlet NSTextField *imageWidthTF;
	IBOutlet NSTextField *imageHeightTF;
	
	IBOutlet OpenGLViewPort * openGLViewport;
	
	IBOutlet NSMatrix *featureOperationType;
	IBOutlet NSButton *useRepetitionButton;
	
	IBOutlet NSMatrix *matchingOperationType;
	IBOutlet NSButton *useMaskButton;
	
	IBOutlet NSMatrix *bundleOperationType;
	
	IBOutlet NSMatrix *visualizationType;
	
	IBOutlet NSTextField * totalNoImagesTF;
	
	IBOutlet NSTextField * currentMatchImageIndexTF;
	IBOutlet NSTextField * neighborImageIndexTF;
	IBOutlet NSTextField * noMatchesTF;
	IBOutlet NSTextField * currentMatchIndexTF;
	
	IBOutlet NSTextField * totalNoTracksTF;
	IBOutlet NSTextField * currentTrackIndexTF;
	IBOutlet NSTextField * noVisibleImagesTF;
	IBOutlet NSTextField * currentTrackImageIndexTF;
	
	IBOutlet NSTextField * totalNoPlanesTF;
	
	IBOutlet NSButton *occluderCheckBox;
	IBOutlet NSMatrix *propagationType;
	
	IBOutlet NSButton *addTextureCheckBox;
	IBOutlet NSButton *repeatTextureButton;
	
	VISUALIZATION_TYPE visType;
	float imgScale;
	
	int neighborImage;
	int currentMatch;
	int totalNoMatches;
	
	int totalNoTracks;
	int totalNoImages;
	int visibleImageNo;
	int currentTrack;
	int currentImage;
	
	int currentPlane;
	int numberOfPlanes;
}

- (void) updateFeatureInfo;
- (void) updateRepetitionInfo;
- (void) updateFeatureMatchInfo;
- (void) updateTrackInfo;
- (void) updateGridStrokeInfo;

- (IBAction)openDocument:(id)sender;

- (IBAction) runFeatureOperation:(id)sender;
- (IBAction) runMatchingOperation:(id)sender;
- (IBAction) runBundleOperation:(id)sender;

- (IBAction)featurePointVisualization:(id)sender;
- (IBAction)repetitionVisualization:(id)sender;
- (IBAction)featureMatchVisualization:(id)sender;
- (IBAction)trackCorrespondenceVisualization:(id)sender;
- (IBAction)gridStrokeVisualization:(id)sender;

- (IBAction)nextNeighbor:(id)sender;
- (IBAction)prevNeighbor:(id)sender;
- (IBAction)nextMatch:(id)sender;
- (IBAction)prevMatch:(id)sender;
- (IBAction)nextTrack:(id)sender;
- (IBAction)prevTrack:(id)sender;
- (IBAction)nextPlane:(id)sender;
- (IBAction)prevPlane:(id)sender;
- (IBAction)nextImage:(id)sender;
- (IBAction)prevImage:(id)sender;

- (IBAction) rectifyImage:(id)sender;
- (IBAction) computeEdges:(id)sender;

- (IBAction)addManualFeaturePoint:(id)sender;
- (IBAction)markMatchUnvalid:(id)sender;
- (IBAction)markMatchValid:(id)sender;

- (IBAction) markTemplate:(id)sender;
- (IBAction) findGridMatches:(id)sender;
- (IBAction) readGridMatches:(id)sender;
- (IBAction) registerImages:(id)sender;
- (IBAction) filterAlignments:(id)sender;
- (IBAction) read3DGridInfo:(id)sender;
- (IBAction) form3DGrid:(id)sender;

- (IBAction) fitPlanesTo3DData:(id)sender;
- (IBAction) rectifyImageWith3DPlane:(id)sender;
- (IBAction) updateRepetitions:(id)sender;

- (IBAction) projectSketchToImages:(id)sender;

- (IBAction) takeScreenshot:(id)sender;
- (IBAction) translateImageLeft:(id)sender;
- (IBAction) translateImageRight:(id)sender;
- (IBAction) translateImageUp:(id)sender;
- (IBAction) translateImageDown:(id)sender;
- (IBAction) scaleImage:(id)sender;

@property (assign) IBOutlet NSWindow *window;

@end
