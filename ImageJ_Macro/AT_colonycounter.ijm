// imageJ macro to count colonies on agar plates

// variable names of the windows
original = getTitle();
background = original + "_bg";

// preliminary formatting on original image
run("8-bit");
run("Invert");

// cropping out the outside of the petri dish as necessary
// dictated by user


// creation of the background which will substracted
run("Duplicate...", "title=["+background+"]");
selectWindow(background);
run("Gaussian Blur...", "sigma=50");

// creating the parsed image
imageCalculator("Subtract create", original, background);
result = substring(original, 0,(lengthOf(original) - 4));
run("Rename...", "title=["+result+"]");


// closing unneeded window
selectWindow(background);
run("Close");

// initial parsing/cleaning of the result
selectWindow(result);
run("Median...", "radius=4");
resetMinAndMax();
run("Enhance Contrast", "saturated=0.35");

// user sets the threshold manually, clicks "ok" when finished
run("Threshold...");
waitForUser("Set threshold, apply, and click OK.");
selectWindow("Threshold");
run("Close");

// splits some of the larger grouped particles (touching colonies)
run("Watershed");

// analyzing
run("Analyze Particles...", "size=0-0.005 circularity=0.40-1.00 show=[Outlines] display clear summarize");

// closing unnecessary windows
selectWindow("Results");
run("Close");
selectWindow("Drawing of " + result);
waitForUser("Close window of colonies counted.");
run("Close");
selectWindow(result);
run("Close");

selectWindow(original);
run("Invert");
