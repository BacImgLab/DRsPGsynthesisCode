// ---------------------------
// Batch ThunderSTORM Analysis Macro
// ---------------------------
// This macro recursively processes folders containing SMT tracking image stacks
// (named like "trackXXX.tif") and performs ThunderSTORM single-molecule localization.
// It saves the localization results as .csv files in each corresponding folder.
// ---------------------------

// Prompt user to choose the top-level folder
inputFolder = getDirectory("ThunderSTORM - Choose the input folder!");

// Get a list of all items (files and subfolders) in the selected directory
subfolders = getFileList(inputFolder);

// Loop through each item in the top-level directory
for (i = 0; i < subfolders.length; i++) {
    subfolder = subfolders[i];

    // Skip system files "." and ".."
    if (subfolder == "." || subfolder == "..") continue;

    // Build full path to the current subfolder
    subfolderPath = inputFolder + subfolder + "/";

    // Get all files inside the subfolder
    images = getFileList(subfolderPath);

    // Loop through all files in the subfolder
    for (j = 0; j < images.length; j++) {
        imageName = images[j];

        // Process only .tif files that start with "track"
        if (startsWith(imageName, "track") && endsWith(imageName, ".tif")) {

            // Full path of the image file
            imagePath = subfolderPath + imageName;

            // Define the output CSV filename by replacing ".tif" with ".csv"
            resultCSVPath = subfolderPath + replace(imageName, ".tif", ".csv");

            // --- Begin ThunderSTORM Processing ---
            // Close all open images and windows to avoid memory issues
            run("Close All");
            wait(500);  // Brief pause to ensure ImageJ is ready

            // Open the SMT image stack for localization
            open(imagePath);
            wait(1000);  // Wait for image to fully load

            // Set up ThunderSTORM camera settings
            // These values depend on your EMCCD/CMOS camera and pixel size
            run("Camera setup", 
                "offset=0 " +
                "quantumefficiency=0.95 " +
                "isemgain=true " +
                "photons2adu=0.59 " +
                "gainem=1 " +
                "pixelsize=110");

            wait(1000);  // Wait for settings to apply

            // Run the core ThunderSTORM analysis pipeline
            run("Run analysis", 
                "filter=[Wavelet filter (B-Spline)] scale=2.0 order=3 " +
                "detector=[Local maximum] connectivity=8-neighbourhood threshold=1.2*std(Wave.F1) " +
                "estimator=[PSF: Gaussian] sigma=1.5 fitradius=3 method=[Least squares] " +
                "full_image_fitting=false mfaenabled=false renderer=[No Renderer]");

            wait(2000);  // Wait for analysis to complete

            // Optional: filter out poor-quality localizations (low intensity, large sigma)
            run("Show results table", 
                "action=filter formula=[sigma < 250 & intensity > 120 & intensity < 700]");

            wait(2000);  // Allow filter to complete

            // Export results to CSV with selected metadata
            run("Export results", 
                "filepath=[" + resultCSVPath + "] " +
                "fileformat=[CSV (comma separated)] " +
                "id=true frame=true sigma=true chi2=true uncertainty_xy=true " +
                "bkgstd=true intensity=true saveprotocol=true offset=true y=true x=true");

            wait(2000);  // Ensure save completes

            // Clean up by closing all windows
            run("Close All");
        }
    }
}
