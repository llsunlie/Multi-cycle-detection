#include "MyMulti-cycle detection.h"

#define FILE_PATH "multi-cycle3z.bmp"
#define TEST_FILE_PATH "testOutput.txt"
#define OUTPUT_FILE_PATH "finalBMP.bmp"

int main(){
	// init bmp
	BMP bmp(FILE_PATH);
	
	// init MyMulticycleDection
	MyMulticycleDetection myMulticycleDetection(bmp);
	// run!!!
	myMulticycleDetection.run();
	// BMP to File
	myMulticycleDetection.finalBMPToFile(bmp, OUTPUT_FILE_PATH);
    printf("done!\n");
	

    return 0;
}