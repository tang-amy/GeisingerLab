// Yunfei Dai
// 01/10/2020
// This macro substractes background for all images open, and concatenate them into hyperstacked TIFF
// To bypass pop-up window: Plugins -> Bio-Formats -> Bio-Formats Plugins Configuration -> Formats -> Choose "Zeiss CZI" -> Check "Windowless"
 

IN_DIR = getDirectory("Choose a Directory");
//The input directory contains multiple folders, each contains just czi images
OUT_DIR = IN_DIR;

list = getFileList(IN_DIR);
for (i = 0; i < list.length; i++) {
                image_dir = IN_DIR + list[i];
                if(File.isDirectory(image_dir))
                	BGredc(image_dir, OUT_DIR);
}

run("Bio-Formats Macro Extensions")
setBatchMode(true);

function BGredc(input, output) {
	images = getFileList(input);
	title = File.getNameWithoutExtension(images[0]);
	outname = output + title + "BGredc.tif" ;
	commandString = "title=" + title + " open " ;
	for (i = 0; i < images.length; i++) {
		fname = input + images[i];
		commandString = commandString + "image"+(i+1) + "=" + images[i] + " ";
        open(fname);
		run("Subtract Background...", "rolling=50 slice");
	}
	run("Concatenate...", commandString);
	saveAs("Tiff", outname);
	close("*");
}


