//
//  BundlerTrackViewerAppDelegate.m
//  BundlerTrackViewer
//
//  Created by Duygu Ceylan on 11/8/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#import "BundlerTrackViewerAppDelegate.h"
#include "BundlerOutputParser.h"
#include "BundlerManager.h"

@implementation BundlerTrackViewerAppDelegate

@synthesize window;

- (void)applicationDidFinishLaunching:(NSNotification *)aNotification {
	// Insert code here to initialize your application 
	bundlerManager = NULL;
	visType = FeaturePoints;
	imgScale = 1.0;
}

- (IBAction)openDocument:(id)sender
{
	// Create the File Open Dialog class.
	NSOpenPanel* openDlg = [NSOpenPanel openPanel];
	
	// Enable the selection of files in the dialog.
	[openDlg setCanChooseFiles:YES];
	
	// Display the dialog.  If the OK button was pressed,
	// process the files.
	if ( [openDlg runModalForDirectory:nil file:nil] == NSOKButton )
	{
		// Get an array containing the full filenames of all
		// files and directories selected.
		const char *fname = [[openDlg filename] UTF8String];
		NSString *name = [[NSString alloc] initWithCString:fname];
		
		NSString *fileWrapperPath;
		const char *mainDirectory;
		
		//compare extension
		if([name	hasSuffix:[[NSString alloc] initWithCString:".png"]] || [name	hasSuffix:[[NSString alloc] initWithCString:".jpg"]])
		{
			totalNoImages = 0;
			mainDirectory = [[[name stringByDeletingLastPathComponent] stringByDeletingLastPathComponent] UTF8String];
			fileWrapperPath = [name stringByDeletingLastPathComponent];
		}
		if([name	hasSuffix:[[NSString alloc] initWithCString:".out"]])
		{
			//images
			totalNoImages = 0;
			fileWrapperPath = [[name stringByDeletingLastPathComponent] stringByDeletingLastPathComponent];
			mainDirectory = [fileWrapperPath UTF8String];
			fileWrapperPath = [fileWrapperPath stringByAppendingString:@"/images"];
		}
		
		bundlerManager = new BundlerManager(string(mainDirectory));
		
		NSLog(@"fileWrapperPath = %@", fileWrapperPath);
			
		// load image sequence
		NSFileWrapper * fileWrapper = [[NSFileWrapper alloc] initWithPath:fileWrapperPath];
			
		[fileWrapper autorelease];
			
		if([fileWrapper isDirectory])
		{
			NSDictionary * fileWrapperContent = [fileWrapper fileWrappers];
				
			NSLog(@"File wrapper is a directory with %d files",[fileWrapperContent count]);
				
			NSArray * listOfAllFilePaths = [fileWrapperContent allKeys];
			NSMutableArray * listOfFilePaths = [[NSMutableArray alloc] init];
				
			for(unsigned int i=0;i<[listOfAllFilePaths count];i++)
			{
				NSString * filePathString = [listOfAllFilePaths objectAtIndex:i];
				NSRange rangePng = [filePathString rangeOfString:@".png"];
				NSRange rangeJpg = [filePathString rangeOfString:@".jpg"];
				if(rangePng.location != NSNotFound || rangeJpg.location != NSNotFound) 
				{
					[listOfFilePaths addObject:[listOfAllFilePaths objectAtIndex:i]];
				}
			}
				
			NSArray * filePathsSorted = [listOfFilePaths sortedArrayUsingSelector:@selector(compare:)];
			for(unsigned int i=0;i<[filePathsSorted count];i++)
			{
				NSString * key = [filePathsSorted objectAtIndex:i];
				NSString * filePath = [NSString stringWithString:fileWrapperPath];
				filePath = [filePath stringByAppendingString:@"/"];
				filePath = [filePath stringByAppendingString:key];
				
				NSLog(@"image filename of index %d:%@",i,filePath);
					
				std::string stdFilePath = std::string([filePath UTF8String]);
				bundlerManager->addNewImagePath(stdFilePath);
				[openGLViewport addNewImagePath:stdFilePath];
				totalNoImages++;
			}
		}
		
		//check if mask images exist
		fileWrapperPath = [[NSString stringWithUTF8String:mainDirectory] stringByAppendingString:@"/masks"];
		NSFileWrapper * maskWrapper = [[NSFileWrapper alloc] initWithPath:fileWrapperPath];
		
		[maskWrapper autorelease];
		
		if([maskWrapper isDirectory])
		{
			NSDictionary * maskWrapperContent = [maskWrapper fileWrappers];
			
			NSLog(@"Mask wrapper is a directory with %d files",[maskWrapperContent count]);
			
			NSArray * listOfAllFilePaths = [maskWrapperContent allKeys];
			NSMutableArray * listOfFilePaths = [[NSMutableArray alloc] init];
			
			for(unsigned int i=0;i<[listOfAllFilePaths count];i++)
			{
				NSString * filePathString = [listOfAllFilePaths objectAtIndex:i];
				NSRange rangePng = [filePathString rangeOfString:@".png"];
				NSRange rangeJpg = [filePathString rangeOfString:@".jpg"];
				if(rangePng.location != NSNotFound || rangeJpg.location != NSNotFound) 
				{
					[listOfFilePaths addObject:[listOfAllFilePaths objectAtIndex:i]];
				}
			}
			
			NSArray * filePathsSorted = [listOfFilePaths sortedArrayUsingSelector:@selector(compare:)];
			for(unsigned int i=0;i<[filePathsSorted count];i++)
			{
				NSString * key = [filePathsSorted objectAtIndex:i];
				NSString * filePath = [NSString stringWithString:fileWrapperPath];
				filePath = [filePath stringByAppendingString:@"/"];
				filePath = [filePath stringByAppendingString:key];
				NSLog(@"mask filename of index %d:%@",i,filePath);
				
				std::string stdFilePath = std::string([filePath UTF8String]);
				bundlerManager->addNewMaskImagePath(stdFilePath);
			}
		}
		
		[totalNoImagesTF setIntValue:totalNoImages];
		
		[openGLViewport setScaleOfImg:imgScale];
		
		if(totalNoImages > 0)
		{
			//[openGLViewport loadAllImages:scale];
			[openGLViewport updateActiveImage:0];
			[openGLViewport display];
		}
		
		if([name	hasSuffix:[[NSString alloc] initWithCString:".out"]])
		{
			bundleParser = new BundlerOutputParser();
			[openGLViewport attachBundlerParser:bundleParser];
			
			int imageWidth = [imageWidthTF intValue];
			int imageHeight = [imageHeightTF intValue];
			bundleParser->readBundlerOutputFile(fname, imageWidth, imageHeight);
			totalNoTracks = bundleParser->getNoOfTrack();
			[totalNoTracksTF setIntValue: totalNoTracks];
			//first track
			currentTrack = 0;
			currentImage = 0;
			[self updateTrackInfo];
		}
	}		
}

-(IBAction) runFeatureOperation:(id)sender
{
	NSButtonCell *selCell = [featureOperationType selectedCell];
	NSString *title = [selCell title];
	
	if([title compare:@"Compute Feature Points"] == NSOrderedSame)
	{
		if(bundlerManager != NULL)
			bundlerManager->ComputeAllFeaturePoints();
		
		if(visType == FeaturePoints)
		{
			[self updateFeatureInfo];
		}
		
	}
	else if([title compare:@"Find Unique Feature Points"] == NSOrderedSame)
	{
		if(bundlerManager != NULL)
			bundlerManager->findUniqueFeaturePoints();
		
		if(visType == FeaturePoints)
		{
			[self updateFeatureInfo];
		}
	}
	else if([title compare:@"Save Feature Points"] == NSOrderedSame)
	{
		if(bundlerManager != NULL)
			bundlerManager->saveFeaturePoints();
	}
	else if([title compare:@"Read Feature Points"] == NSOrderedSame)
	{
		if(bundlerManager != NULL)
		{
			bool useRepetition = [useRepetitionButton state];
			
			if(useRepetition)
				bundlerManager->ReadAllFeaturePointsWithRepetition();
			else
				bundlerManager->ReadAllFeaturePoints();
			
			if(visType == FeaturePoints)
			{
				[self updateFeatureInfo];
			}
		}
	}
	
}

-(IBAction) runMatchingOperation:(id)sender
{
	NSButtonCell *selCell = [matchingOperationType selectedCell];
	NSString *title = [selCell title];
	
	if([title compare:@"Match Feature Points"] == NSOrderedSame)
	{
		if(bundlerManager != NULL)
		{
			bool useMask = [useMaskButton state];
			bundlerManager->matchFeaturePoints(useMask);
		}
	}
	else if([title compare:@"Build Manual Matches"] == NSOrderedSame)
	{
		if(bundlerManager != NULL)
			bundlerManager->buildManualCorrespondences();
	}
	else if([title compare:@"Save Feature Matches"] == NSOrderedSame)
	{
		if(bundlerManager != NULL)
			bundlerManager->writeMatches();
	}
	else if([title compare:@"Read Feature Matches"] == NSOrderedSame)
	{
		if(bundlerManager != NULL)
			bundlerManager->readMatches();
	}
	
}

-(IBAction) runBundleOperation:(id)sender
{
	NSButtonCell *selCell = [bundleOperationType selectedCell];
	NSString *title = [selCell title];
	
	if([title compare:@"Bundle Adjustment"] == NSOrderedSame)
	{
		if(bundlerManager != NULL)
		{
			bundlerManager->runBundleAdjustment();
		}
	}
	else if([title compare:@"Save Calibration"] == NSOrderedSame)
	{
		if(bundlerManager != NULL)
		{
			bundlerManager->writeCalibrationResults();
		}
	}
	else if([title compare:@"Read Calibration"] == NSOrderedSame)
	{
		if(bundlerManager != NULL)
		{
			bundlerManager->readCalibrationResults();
		}
	}
}

- (IBAction)featurePointVisualization:(id)sender
{
	visType = FeaturePoints;
	[openGLViewport updateActiveImage:currentImage];
	[self updateFeatureInfo];
}

- (IBAction) featureMatchVisualization:(id)sender
{
	visType = FeatureMatches;
	currentImage = 0;
	neighborImage = 1;
	currentMatch = 0;
	[openGLViewport updateActiveImages:currentImage linkedTo:neighborImage];
	
	[self updateFeatureMatchInfo];
}

- (IBAction) repetitionVisualization:(id)sender
{
	visType = Repetitions;
	[openGLViewport updateActiveImage:currentImage];
	[self updateRepetitionInfo];
}

- (IBAction) trackCorrespondenceVisualization:(id)sender
{
	visType = Tracks;
	if(bundlerManager != NULL)
	{
		totalNoTracks = bundlerManager->getNoOfTrack();
		[totalNoTracksTF setIntValue: totalNoTracks];
		//first track
		currentTrack = 0;
		currentImage = 0;
		[self updateTrackInfo];
	}
}

- (IBAction) gridStrokeVisualization:(id)sender
{
	visType = GridStrokes;
	if(bundlerManager != NULL)
	{
		currentImage = 0;
		[openGLViewport updateActiveImage:currentImage];
		[self updateGridStrokeInfo];
	}
}

-(void) updateFeatureInfo
{
	keypt_t *featurePointInfo;
	int noOfFeatures;
	bundlerManager->getFeaturePointInfo(currentImage, featurePointInfo, noOfFeatures);
	[openGLViewport setFeaturePointInfo:featurePointInfo withFeatureNo:noOfFeatures];
	[openGLViewport updateActiveImage:currentImage];
	int noOfProjPoints;
	float *projections = bundlerManager->getProjectedPoints(noOfProjPoints, currentImage);
	Color *projColors = bundlerManager->getProjectedPointColors(noOfProjPoints, currentImage);
	
	[openGLViewport setProjectedPointInfo:projections withColor:projColors withPointNo:noOfProjPoints];
	[openGLViewport display];
	
}

-(void) updateFeatureMatchInfo
{
	totalNoMatches = bundlerManager->getNoMatches(currentImage, neighborImage);
	[currentMatchImageIndexTF setIntValue:currentImage];
	[neighborImageIndexTF setIntValue:neighborImage];
	[noMatchesTF setIntValue:totalNoMatches];
	
	Vec2f pt1, pt2;
	Vec2f epipolarStart, epipolarEnd;
	bool valid;
	bundlerManager->getMatchingPoints(currentImage, neighborImage, currentMatch, pt1, pt2, epipolarStart, epipolarEnd, valid);
	pt1 /= imgScale; pt2 /= imgScale;
	[currentMatchIndexTF setIntValue:currentMatch];
	if(epipolarStart != Vec2f(0.0, 0.0) && epipolarEnd != Vec2f(0.0, 0.0))
	{
		[openGLViewport updateActiveCorrespondences:pt1 linkedTo:pt2 withEpipolarStart:epipolarStart withEpipolarEnd:epipolarEnd];
	}
	else
	{
		[openGLViewport updateActiveCorrespondences:pt1 linkedTo:pt2 isValid:valid];
	}
	
	/*vector<Vec2f> pt1;
	vector<Vec2f> pt2;
	vector<bool> valid;
	
	bundlerManager->getAllMatchingPoints(currentImage, neighborImage, pt1, pt2, valid);
	[openGLViewport updateCorrespondences:pt1 linkedTo:pt2 validity:valid];
	*/
	
	[openGLViewport display];
}

-(void) updateRepetitionInfo
{
	[openGLViewport updateActiveImage:currentImage];
	//int noPlanes;
	//Img *r = bundlerManager->RectifyImage(currentImage, noPlanes);
	//if(r != NULL)
	//{
	//	[openGLViewport updateRectifiedImage:r];
	//}
	vector<vector<Vec2f> > repetitions = bundlerManager->getRepetitions(currentImage);
	vector<Color> colors = bundlerManager->getRepetitionColors(currentImage);
	vector<vector<Vec2f> > refinedRepetitions = bundlerManager->getRefinedRepetitions(currentImage);
	
	[openGLViewport setRepetitions:repetitions withColors:colors];
	[openGLViewport setRefinedRepetitions:refinedRepetitions];
	
	[openGLViewport display];
	
}

-(void) updateTrackInfo
{
	visibleImageNo = bundleParser->getNoOfCorrespondencesInTrack(currentTrack);
	//visibleImageNo = bundlerManager->getNoOfCorrespondencesInTrack(currentTrack);
	if(visibleImageNo == 0)
		return;
	
	[currentTrackIndexTF setIntValue:currentTrack];
	[noVisibleImagesTF setIntValue:visibleImageNo];
	
	Correspondence c = bundleParser->getCorrespondenceInTrack(currentTrack, currentImage);
	//Correspondence c = bundlerManager->getCorrespondenceInTrack(currentTrack, currentImage);
	[currentTrackImageIndexTF setIntValue:c.imageIndex];
	
	int linkedIndex = currentImage+1;
	if(linkedIndex == visibleImageNo)
		linkedIndex = 0;
	Correspondence cCurrent = bundleParser->getCorrespondenceInTrack(currentTrack, currentImage);
	Correspondence cLinked = bundleParser->getCorrespondenceInTrack(currentTrack, linkedIndex);
	//Correspondence cLinked = bundlerManager->getCorrespondenceInTrack(currentTrack, linkedIndex);
	//cCurrent.position /= imgScale;
	//cLinked.position /= imgScale;
	
	[openGLViewport updateActiveImages:c.imageIndex linkedTo:cLinked.imageIndex];
	[openGLViewport updateActiveCorrespondences:c.position linkedTo:cLinked.position];
	[openGLViewport display];
	//[openGLViewport takeScreenshot];
}

-(void)updateGridStrokeInfo
{
	//vector<vector<Vec2f> > strokes = bundlerManager->getGridStrokes(currentImage);
	//[openGLViewport setGridStrokes:strokes];
	//[openGLViewport display];
}

-(IBAction) nextTrack:(id)sender
{
	currentTrack++;
	if(currentTrack == totalNoTracks)
		currentTrack = 0;
	currentImage = 0;
	[self updateTrackInfo];
}

-(IBAction) prevTrack:(id)sender
{
	currentTrack--;
	if(currentTrack == -1)
		currentTrack = totalNoTracks-1;
	currentImage = 0;
	[self updateTrackInfo];
}

-(IBAction) nextMatch:(id)sender
{
	currentMatch++;
	if(currentMatch >= totalNoMatches)
		currentMatch = 0;
	[self updateFeatureMatchInfo];
}

-(IBAction) prevMatch:(id)sender
{
	currentMatch--;
	if(currentMatch == -1)
		currentMatch = totalNoMatches-1;
	[self updateFeatureMatchInfo];
}

- (IBAction)nextNeighbor:(id)sender
{
	neighborImage++;
	
	if(neighborImage == currentImage)
		neighborImage = currentImage+1;
	
	if(neighborImage == totalNoImages)
		neighborImage = 0;
	
	currentMatch = 0;
	[openGLViewport updateActiveImages:currentImage linkedTo:neighborImage];
	[self updateFeatureMatchInfo];
}

- (IBAction)prevNeighbor:(id)sender
{
	neighborImage--;
	
	if(neighborImage == currentImage)
		neighborImage = currentImage-1;
	if(neighborImage == -1)
		neighborImage = totalNoImages - 1;
	
	currentMatch = 0;
	[openGLViewport updateActiveImages:currentImage linkedTo:neighborImage];
	[self updateFeatureMatchInfo];
}

-(IBAction) nextImage:(id)sender
{
	currentImage++;
	
	if(visType == FeaturePoints)
	{
		if(currentImage == totalNoImages)
			currentImage = 0;
		[self updateFeatureInfo];
	}
	else if(visType == Repetitions)
	{
		if(currentImage == totalNoImages)
			currentImage = 0;
		[self updateRepetitionInfo];
	}
	else if(visType == FeatureMatches)
	{
		if(currentImage == totalNoImages)
			currentImage = 0;
		
		if(currentImage == 0)
			neighborImage = 1;
		else
			neighborImage = 0;
		
		currentMatch = 0;
		[openGLViewport updateActiveImages:currentImage linkedTo:neighborImage];
		[self updateFeatureMatchInfo];
	}
	else if(visType == Tracks)
	{
		if(currentImage == visibleImageNo)
			currentImage = 0;
		[self updateTrackInfo];
	}
	else if(visType == GridStrokes)
	{
		if(currentImage == totalNoImages)
			currentImage = 0;
		[openGLViewport updateActiveImage:currentImage];
		[self updateGridStrokeInfo];
	}
}

-(IBAction) prevImage:(id)sender
{
	currentImage--;
	
	if(visType == FeaturePoints)
	{
		if(currentImage == -1)
			currentImage = totalNoImages - 1;
		[self updateFeatureInfo];
	}
	else if(visType == Repetitions)
	{
		if(currentImage == -1)
			currentImage = totalNoImages - 1;
		[self updateRepetitionInfo];
	}
	else if(visType == FeatureMatches)
	{
		if(currentImage == -1)
			currentImage = totalNoImages - 1;
		if(currentImage == 0)
			neighborImage = 1;
		else
			neighborImage = 0;
		
		currentMatch = 0;
		[openGLViewport updateActiveImages:currentImage linkedTo:neighborImage];
		[self updateFeatureMatchInfo];
	}
	else if(visType == Tracks)
	{
		if(currentImage == -1)
			currentImage = visibleImageNo-1;
		[self updateTrackInfo];
	}
	else if(visType == GridStrokes)
	{
		if(currentImage == -1)
			currentImage = visibleImageNo-1;
		[openGLViewport updateActiveImage:currentImage];
		[self updateGridStrokeInfo];
	}
}

-(IBAction) rectifyImage:(id)sender
{
	if(bundlerManager != NULL)
	{
		Img *r = bundlerManager->RectifyImage(currentImage, numberOfPlanes);
		currentPlane = 0;
		if(r != NULL)
		{
			[totalNoPlanesTF setIntValue:numberOfPlanes];
			[openGLViewport updateRectifiedImage:r];
			[openGLViewport display];
		}
	}
}

-(IBAction) prevPlane:(id)sender
{
	if(numberOfPlanes == 0)
		return;
	currentPlane -= 1;
	if(currentPlane < 0)
	{
		currentPlane = numberOfPlanes-1;
	}
	Img *r = bundlerManager->getRectifiedImage(currentImage, currentPlane);
	[openGLViewport updateRectifiedImage:r];
	[openGLViewport display];
}

-(IBAction) nextPlane:(id)sender
{
	if(numberOfPlanes == 0)
		return;
	currentPlane += 1;
	if(currentPlane == numberOfPlanes)
	{
		currentPlane = 0;
	}
	Img *r = bundlerManager->getRectifiedImage(currentImage, currentPlane);
	[openGLViewport updateRectifiedImage:r];
	[openGLViewport display];
}

- (IBAction) computeEdges:(id)sender
{
	if(bundlerManager != NULL)
	{
		bundlerManager->computeEdges();
	}
}

- (IBAction) markTemplate:(id)sender
{
	[openGLViewport setDrawLine:true];
	[occluderCheckBox setState:NSOffState];
}

- (IBAction) findGridMatches:(id)sender
{
	if(bundlerManager != NULL)
	{
		Vec2f topLeft;
		int w, h;
		
		[openGLViewport getTemplateCoord:topLeft width:w height:h];
		bundlerManager->findMatches(currentImage, currentPlane, topLeft, w, h);
	}
}

- (IBAction) readGridMatches:(id)sender
{
	if(bundlerManager != NULL)
	{
		bundlerManager->readAllGridMatches();
	}
}

- (IBAction) registerImages:(id)sender
{
	if(bundlerManager != NULL)
	{
		bundlerManager->registerImages();
	}
}

- (IBAction) filterAlignments:(id)sender
{
	if(bundlerManager != NULL)
	{
		bundlerManager->convertAlignmentsToMatches();
	}
}

- (IBAction) read3DGridInfo: (id)sender
{
	if(bundlerManager != NULL)
	{
		bundlerManager->read3DGrids();
	}
}

- (IBAction) form3DGrid: (id)sender
{
	if(bundlerManager != NULL)
	{
		bool success = bundlerManager->form3DGrids();
		//if(success)
		//	[openGLViewport setGrids3DComputed:true];
	}
}

- (IBAction) fitPlanesTo3DData:(id)sender
{
	if(bundlerManager != NULL)
	{
		bundlerManager->fitPlanesTo3DData();
	}
}

- (IBAction) rectifyImageWith3DPlane:(id)sender
{
	if(bundlerManager != NULL)
	{
		bundlerManager->rectifyWith3DPlane();
	}
}

- (IBAction) updateRepetitions:(id)sender
{
	if(bundlerManager != NULL)
	{
		bundlerManager->updateGridCells();
	}
}

- (IBAction)addManualFeaturePoint:(id)sender
{
	if(bundlerManager != NULL)
	{
		Vec2f c;
		[openGLViewport getCurrentCoord:c];
		bundlerManager->addFeaturePoint(currentImage, c);
		if(visType == FeaturePoints)
			[self updateFeatureInfo];
	}
}

- (IBAction)markMatchUnvalid:(id)sender
{
	if(bundlerManager != NULL)
	{
		bundlerManager->setMatchValidity(currentImage, neighborImage, currentMatch, false);
	}
}

- (IBAction)markMatchValid:(id)sender
{
	if(bundlerManager != NULL)
	{
		bundlerManager->setMatchValidity(currentImage, neighborImage, currentMatch, true);
	}	
}

- (IBAction) projectSketchToImages:(id)sender
{
	/*if(bundlerManager != NULL)
	{
		vector<Vec2f> endpoints;
		[openGLViewport getStrokes:endpoints];
		
		if([occluderCheckBox state] == NSOnState)
		{
			NSButtonCell *selCell = [propagationType selectedCell];
			NSString *title = [selCell title];
			
			int propationType;
			
			if([title compare:@"Dense Reconstruction"] == NSOrderedSame)
			{
				propationType = 0;
				bundlerManager->projectOccluderSketchToImages(currentImage, endpoints, propationType);
			}
			else if([title compare:@"Depth Estimation"] == NSOrderedSame)
			{
				propationType = 1;
				bundlerManager->projectOccluderSketchToImages(currentImage, endpoints, propationType);
			}
			else if([title compare:@"With Repetitions"] == NSOrderedSame)
			{
				bundlerManager->removeOccluderUsingRepetitions(currentImage, endpoints);
			}
			
			[visualizationType setSelectionFrom:4 to:4 anchor:4 highlight:YES];
			visType = GridStrokes;
			[self updateGridStrokeInfo];
		}
		else if([addTextureCheckBox state] == NSOnState)
		{
			if([repeatTextureButton state] == NSOnState)
			{
				bundlerManager->projectGridSketchToImages(currentImage, endpoints);
			}
			else
			{
				bundlerManager->projectTextureSketchToImages(currentImage, endpoints);
			}
			[visualizationType setSelectionFrom:4 to:4 anchor:4 highlight:YES];
			visType = GridStrokes;
			[self updateGridStrokeInfo];
		}
	}*/
}

-(IBAction) takeScreenshot:(id)sender
{
	[openGLViewport takeScreenshot];
	[openGLViewport display];
}

- (IBAction)translateImageLeft:(id)sender
{
	[openGLViewport translateImageLeft];
	[openGLViewport display];
}
- (IBAction) translateImageRight:(id)sender
{
	[openGLViewport translateImageRight];
	[openGLViewport display];
}
- (IBAction) translateImageUp:(id)sender
{
	[openGLViewport translateImageUp];
	[openGLViewport display];
}
- (IBAction) translateImageDown:(id)sender
{
	[openGLViewport translateImageDown];
	[openGLViewport display];
}
- (IBAction) scaleImage:(id)sender
{
	[openGLViewport scaleImage];
	[openGLViewport display];
}

@end
