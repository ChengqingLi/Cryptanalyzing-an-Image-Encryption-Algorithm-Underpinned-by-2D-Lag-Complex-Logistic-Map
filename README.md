# Cryptanalyzing an Image Encryption Algorithm Underpinned by 2D Lag-Complex-Logistic-Map

This code requires OpenCV basic environment. All the C++ code is contained in the `Main` function(include the plain-images generation process), and the `image` folder is responsible for storing the images.

- Main.cpp

`a, b, SumR, SumG, SumB` correspond to the control parameter and key. In the main function, you can read the image stored in `image` folder (not larger than 3\*256\*256), the code will automatically run the attack process, and the cipher-image `Tempchiper.jpg` and attack result `Tempattack.jpg` during this process will be stored in the `image` folder.
