macro "bersenThtest" {
    // Get the title of the original (active) image
    original = getTitle();
    
    // Set parameter ranges for Bernsen thresholding
    minRadius = 10;    // Minimum radius for local thresholding window
    maxRadius = 30;    // Maximum radius
    stepRadius = 5;    // Step size for radius
    
    minParam1 = 15;    // Minimum parameter_1 value (contrast threshold)
    maxParam1 = 15;    // Maximum parameter_1 value (fixed here)
    stepParam1 = 2;    // Step size for parameter_1
    
    param2 = 0;        // Fixed parameter_2 (required by command, unused by Bernsen)
    
    // Loop through radius values
    for (radius = minRadius; radius <= maxRadius; radius += stepRadius) {
        // Loop through parameter_1 values
        for (param1 = minParam1; param1 <= maxParam1; param1 += stepParam1) {
            // Select the original image window
            selectWindow(original);
            
            // Duplicate the original image and assign a descriptive title including parameters
            run("Duplicate...", "title=Radius" + radius + "_Param1_" + param1);
            
            // Apply Bernsen local thresholding on the duplicated image
            run("Auto Local Threshold", "method=Bernsen radius=" + radius + 
                " parameter_1=" + param1 + " parameter_2=" + param2 + " white");
        }
    }
}
