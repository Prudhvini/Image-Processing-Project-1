/*
 B490/B659 Project 1 Skeleton Code    (1/2015)
 
 Be sure to read over the project document and this code (should you choose to use it) before 
 starting to program. 
 

 Compiling:
 A simple console command 'make' will excute the compilation instructions
 found in the Makefile bundled with this code. It will result in an executable
 named p1.

 Running:
 The executable p1 takes commands of the form:
 ./p1 problem_ID input_File ouput_File additional_Arguments
 
 Some examples:

 ./p1 2.1 input.png out.png 
 
 This runs the 'averageGrayscale' function defined in this file and described
 in the project doc Part 2, problem 1 on the input. The output is saved as out.png. 
 
 ----

 ./p1 4.2b input.jpg output.png 0.5

 This runs the Gaussian filter function (gaussianFilter) with sigma = 0.5 on input.jpg.
 This is problem 2b from Part 4 in the project documentation. 

 */

//Link to the header file
#include "CImg.h"
#include <ctime>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <math.h>
#define TWO_PI 6.2831853071795864769252866
//Use the cimg namespace to access the functions easily
using namespace cimg_library;
using namespace std;

//Part 2 - Basic image operations
CImg<double> averageGrayscale(CImg<double> input);
CImg<double> simpleBW(CImg<double> input);
CImg<double> advancedBW(CImg<double> input);

//Part 3 - Adding noise
CImg<double> uniformNoise(CImg<double> input);
CImg<double> gaussianNoise(CImg<double> input, double sigma);
CImg<double> saltAndPepperNoise(CImg<double> input);

//Part 4 - Filtering
CImg<double> filter(CImg<double> input, CImg<double> filter);
CImg<double> meanFilter(CImg<double> input, int filterSize);
CImg<double> gaussianFilter(CImg<double> input, double sigma);
CImg<double> medianFilter(CImg<double> input, int size);

//Part 5 - Quantitative Analysis
double meanSquareError(CImg<double> input, CImg<double> output);

//Part 6 - Speedup
CImg<double> separableKCGaussian(CImg<double> input);
CImg<double> separableKernelConvolution(CImg<double> input, int size);
CImg<double> dynamicProgrammingConvolution(CImg<double> input, int size);
CImg<double> cascadedGuassianFilter(CImg<double> input, int size);

int main(int argc, char **argv)
{
	
	if (argc < 4)
	{
		cout << "Insufficent number of arguments. Please see documentation" << endl;
		cout << "p1 problemID inputfile outputfile" << endl;
		return -1;
	}
	
	char* inputFile = argv[2];
	char* outputFile = argv[3];
	cout << "In: " << inputFile << "  Out: " << outputFile << endl;
	
	CImg<double> input(inputFile);
	
	CImg<double> output;
	
	if (!strcmp(argv[1], "2.1"))
	{
		cout << "# Problem 2.1 - Average Grayscale" << endl;
		if (input.spectrum() != 3)
		{
			cout << "INPUT ERROR: Input image is not a color image!" << endl;
			return -1;
		}
		output = averageGrayscale(input);
		output.save(outputFile);
	}
	else if (!strcmp(argv[1], "2.2a"))
	{
		cout << "# Problem 2.1a - Simple Threshold Black and White" << endl;
		if (input.spectrum() != 1)
		{
			cout << "INPUT ERROR: Input image is not a grayscale image!" << endl;
			return -1;
		}
		output = simpleBW(input);
		output.save(outputFile);
	}
	else if (!strcmp(argv[1], "2.2b"))
	{
		cout << "# Problem 2.2b - Advanced Threshold Black and White" << endl;
		if (input.spectrum() != 1)
		{
			cout << "INPUT ERROR: Input image is not a grayscale image!" << endl;
			return -1;
		}
		output = advancedBW(input);
		output.save(outputFile);
	}
	else if (!strcmp(argv[1], "3.1"))
	{
		cout << "# Problem 3.1 - Uniform Noise" << endl;
		if (input.spectrum() != 1)
		{
			cout << "INPUT ERROR: Input image is not a grayscale image!" << endl;
			return -1;
		}
		output = uniformNoise(input);
		output.save(outputFile);
	}
	else if (!strcmp(argv[1], "3.2"))
	{
		cout << "# Problem 3.2 - Gaussian Noise" << endl;
		if (input.spectrum() != 1)
		{
			cout << "INPUT ERROR: Input image is not a grayscale image!" << endl;
			return -1;
		}
		if (argc != 5)
		{
			cout << "INPUT ERROR: Provide sigma as additional argument!" << endl;
			return -1;
		}
		output = gaussianNoise(input, atof(argv[4]));
		output.save(outputFile);
	}
	else if (!strcmp(argv[1], "3.3"))
	{
		cout << "# Problem 3.3 - Salt & Pepper Noise" << endl;
		if (input.spectrum() != 1)
		{
			cout << "INPUT ERROR: Input image is not a grayscale image!" << endl;
			return -1;
		}
		output = saltAndPepperNoise(input);
		output.save(outputFile);
	}
	else if (!strcmp(argv[1], "4.2a"))
	{
		cout << "# Problem 4.2a - Mean Filter Noise" << endl;
		if (input.spectrum() != 1)
		{
			cout << "INPUT ERROR: Input image is not a grayscale image!" << endl;
			return -1;
		}
		if (argc != 5)
		{
			cout << "INPUT ERROR: Provide filter size as additional argument!" << endl;
			return -1;
		}
		std::clock_t start = std::clock();
		output = meanFilter(input, atoi(argv[4]));
		output.save(outputFile);
		std::cout << "Time Elapsed: " << (std::clock() - start) / (double) (CLOCKS_PER_SEC / 1000) << " ms" << std::endl;
	}
	else if (!strcmp(argv[1], "4.2b"))
	{
		cout << "# Problem 4.2b - Gaussian Filter Noise" << endl;
		if (input.spectrum() != 1)
		{
			cout << "INPUT ERROR: Input image is not a grayscale image!" << endl;
			return -1;
		}
		if (argc != 5)
		{
			cout << "INPUT ERROR: Provide sigma as additional argument!" << endl;
			return -1;
		}
		std::clock_t start = std::clock();
		output = gaussianFilter(input, atof(argv[4]));
		output.save(outputFile);
		std::cout << "Time Elapsed: " << (std::clock() - start) / (double) (CLOCKS_PER_SEC / 1000) << " ms" << std::endl;
	}
	else if (!strcmp(argv[1], "4.3"))
	{
		cout << "# Problem 4.3 - Median Noise" << endl;
		if (input.spectrum() != 1)
		{
			cout << "INPUT ERROR: Input image is not a grayscale image!" << endl;
			return -1;
		}
		if (argc != 5)
		{
			cout << "INPUT ERROR: Provide filter size as additional argument!" << endl;
			return -1;
		}
		output = medianFilter(input, atoi(argv[4]));
		output.save(outputFile);
	}
	else if (!strcmp(argv[1], "5"))
	{
		CImg<double> input2(outputFile);
		cout << "# Problem 5 - Noise Removal Analysis" << endl;
	
		if (input.spectrum() != 1)
		{
			cout << "INPUT ERROR: Input image is not a grayscale image!" << endl;
			return -1;
		}
		double error = meanSquareError(input,input2);
		cout<<"\n Mean Square Error = "<<error;
	}
	else if (!strcmp(argv[1], "6.1a"))
	{
		cout << "# Problem 6.1a - Separable Mean Filter Kernel Convolutions" << endl;
		
		std::clock_t start = std::clock();
		
		if (input.spectrum() != 1)
		{
			cout << "INPUT ERROR: Input image is not a grayscale image!" << endl;
			return -1;
		}
		if (argc != 5)
		{
			cout << "INPUT ERROR: Provide filter size as additional argument!" << endl;
			return -1;
		}
		
		output = separableKernelConvolution(input, atoi(argv[4]));
		output.save(outputFile);
	
		std::cout << "Time Elapsed: " << (std::clock() - start) / (double) (CLOCKS_PER_SEC / 1000) << " ms" << std::endl;
		
	}
	else if (!strcmp(argv[1], "6.1b"))
	{
		cout << "# Problem 6.1b - Separable Guassian Kernel Convolutions" << endl;
		std::clock_t start = std::clock();
		output = separableKCGaussian(input);
		output.save(outputFile);
	
		std::cout << "Time Elapsed: " << (std::clock() - start) / (double) (CLOCKS_PER_SEC / 1000) << " ms" << std::endl;
	}
	else if (!strcmp(argv[1], "6.2"))
	{
		cout << "# Problem 6.2 - Dynamic Box Filter" << endl;
		
		std::clock_t start = std::clock();
		if (input.spectrum() != 1)
		{
			cout << "INPUT ERROR: Input image is not a grayscale image!" << endl;
			return -1;
		}
		if (argc != 5)
		{
			cout << "INPUT ERROR: Provide filter size as additional argument!" << endl;
			return -1;
		}
		
		output = dynamicProgrammingConvolution(input, atoi(argv[4]));
		output.save(outputFile);		
		std::cout << "Time Elapsed: " << (std::clock() - start) / (double) (CLOCKS_PER_SEC / 1000) << " ms" << std::endl;
		
	}
	else if (!strcmp(argv[1], "6.3"))
	{
		cout << "# Problem 6.3 - Fast Gaussian Smoothing" << endl;
		
		std::clock_t start = std::clock();
		if (input.spectrum() != 1)
		{
			cout << "INPUT ERROR: Input image is not a grayscale image!" << endl;
			return -1;
		}
		if (argc != 5)
		{
			cout << "INPUT ERROR: Provide filter size as additional argument!" << endl;
			return -1;
		}
		output = cascadedGuassianFilter(input, atoi(argv[4]));
		output.save(outputFile);		
		
		std::cout << "Time Elapsed: " << (std::clock() - start) / (double) (CLOCKS_PER_SEC / 1000) << " ms" << std::endl;
	}
	else
	{
		cout << "Unknown input command" << endl;
	}
	
	return 0;
}

//Part 2 - Basic image operations
CImg<double> averageGrayscale(CImg<double> input)
{
	//Creates a new grayscale image with same size as the input initialized to all 0s (black)
	CImg<double> output(input.width(), input.height(), 1, 1);
	
	for (int i = 0; i < input.height(); i++)
	{
		for (int j = 0; j < input.width(); j++)
		{
			const double sum = (0.3 * input(j, i, 0, 0)) + (0.6 * input(j, i, 0, 1)) + (0.1 * input(j, i, 0, 2));
			output(j, i, 0, 0) = sum;
		}
	}
	return output;
}

CImg<double> averageGrayscale1(CImg<double> input)
{
	double sigma = 20.0;
	int filterSize = 60;
	//Creates a new grayscale image with same size as the input initialized to all 0s (black)
	CImg<double> H(60, 60, 1, 1);
	
	for (int i = 0; i < filterSize; i++)
	{
		for (int j = 0; j < filterSize; j++)
		{
			H(j - ((filterSize + 1) / 2), i - ((filterSize + 1) / 2)) = (exp((-1 * (i * i + j * j)) / (2 * sigma * sigma))) / (2 * M_PI * sigma * sigma);
		}
	}
	
	return H;
}

CImg<double> simpleBW(CImg<double> input)
{
	//Creates a new grayscale image with same size as the input initialized to all 0s (black)
	CImg<double> output(input, "xyzc", 0);
	
	double threshold = 127;
	for (int i = 0; i < input.height(); i++)
	{
		for (int j = 0; j < input.width(); j++)
		{
			if (input(j, i) >= threshold)
				output(j, i) = 255;
			else
				output(j, i) = 0;
		}
	}
	return output;
}

CImg<double> advancedBW(CImg<double> input)
{
	//Creates a new grayscale image with same size as the input initialized to all 0s (black)
	CImg<double> output(input.width(), input.height(), 1, 1);
	
	double threshold = 127;
	double difference = 0;
	double error = 0;
	for (int i = 1; i < input.width() - 1; i++)
	{
		for (int j = 1; j < input.height() - 1; j++)
		{
			if (input(i, j) >= threshold)
			{
				output(i, j) = 255;
				error = input(i, j) - output(i, j);
				input(i, j + 1) += 0.4 * error;
				input(i + 1, j) += 0.3 * error;
				input(i + 1, j - 1) += 0.2 * error;
				input(i + 1, j + 1) += 0.1 * error;
			}
			else
			{
				output(i, j) = 0;
				error = input(i, j) - output(i, j);
				input(i, j + 1) += 0.4 * error;
				input(i + 1, j) += 0.3 * error;
				input(i + 1, j - 1) += 0.2 * error;
				input(i + 1, j + 1) += 0.1 * error;
			}
		}
	}
	
	return output;
}

//Part 3 - Adding noise
CImg<double> uniformNoise(CImg<double> input)
{
	
	//Creates a new image with same size as the input initialized to all 0s (black)
	CImg<double> output(input.width(), input.height(), 1, 1);
	int min = -20;
	int max = 20;
	int randNum;
	for (int i = 0; i < input.width(); i++)
	{
		for (int j = 0; j < input.height(); j++)
		{
			randNum = rand() % (max - min + 1) + min;
			output(i, j) = input(i, j) + randNum;
			if (output(i, j) < 0)
				output(i, j) = 0;
			if (output(i, j) > 255)
				output(i, j) = 255;
		}
	}
	
	return output;
}

double getGaussianNoise(double sigma)
{
    int MAX_VALUE = 3.0 * sigma; //generally 99.5% of gaussian is within 3 standard deviations
    int MIN_VALUE = -1.0 * 3.0 * sigma; //generally 99.5% of gaussian is within 3 standard deviations
    double sigmasq = sigma * sigma;
    double root2pi = sqrt(2 * M_PI);
    double z = rand() % (MAX_VALUE - MIN_VALUE + 1) + MIN_VALUE;
    
    int i = rand() % 2;
    
    if (i == 0)
    {
        return exp((z * z) / (-2 * sigmasq)) / (sigma * root2pi);
    }
    else
    {
        return exp((z * z) / (-2 * sigmasq)) / (-1 * sigma * root2pi);
    }
}

CImg<double> gaussianNoise(CImg<double> input, double sigma)
{
    //Creates a new image with same size as the input initialized to all 0s (black)
    
    CImg<double> output(input.width(), input.height(), 1, 1);
    
    int count = 0;
    for (int i = 0; i < input.width(); i++)
    {
        for (int j = 0; j < input.height(); j++)
        {
            double noise = getGaussianNoise(sigma);
            output(i, j) = input(i, j) + noise * 100;
            if (output(i, j) < 0)
            {
                output(i, j) = 0;
            }
            if (output(i, j) > 255)
            {
                output(i, j) = 255;
            }
        }
    }
    return output;
}

CImg<double> saltAndPepperNoise(CImg<double> input)
{
	
	//Creates a new image with same size as the input initialized to all 0s (black)
	CImg<double> output(input.width(), input.height(), 1, 1);
	double randNum;
	int min = -20;
	int max = 20;
	for (int i = 0; i < input.width(); i++)
	{
		for (int j = 0; j < input.height(); j++)
		{
			do
			{
				randNum = rand() % (max - min + 1) + min;
			} while (randNum == 0);
			
			if (randNum == -20)
				output(i, j) = 0;
			else if (randNum == 20)
				output(i, j) = 255;
			else
				output(i, j) = input(i, j);
		}
	}

	return output;
}

//Part 4 - Filtering

CImg<double> fixBorder(CImg<double> output, int filterSize)
{
	for(int j=filterSize/2;j>0;j--)
	{
		for(int i=0;i<output.height();i++)
			output(j-1,i) = output(j,i);

		for(int i=0;i<output.width();i++)
			output(i,j-1) = output(i,j);
		
		for(int i=0;i<output.height();i++)
			output(output.width()-j,i) = output(output.width()-(j+1),i);
	
		for(int i=0;i<output.width();i++)
			output(i,output.height()-j) = output(i,output.height()-(j+1));
	}
	return output;
}


void printKernel(CImg<double> filter)
{
    for (int i = 0; i < filter.width(); i++)
    {
        for (int j = 0; j < filter.height(); j++)
        {
            cout << filter(i, j) << "\t";
        }
        cout << endl;
    }
    cout << endl;
}

CImg<double> flipFilter(CImg<double> filter)
{
    CImg<double> x(filter, "xyzc", 0);
    
    //flip horizontally
    for (int i = 0; i < filter.width(); i++)
    {
        for (int j = 0; j < filter.height(); j++)
        {
            x(filter.width() - i - 1, j) = filter(i, j);
        }
    }

    //create a copy of horizontally flipped x
    CImg<double> y(x, "xyzc", 0);
    
    for (int i = 0; i < filter.width(); i++)
    {
        for (int j = 0; j < filter.height(); j++)
        {
            y(i, j) = x(i, j);
        }
    }
    
    //flip vertically
    for (int i = 0; i < filter.height(); i++)
    {
        for (int j = 0; j < filter.width(); j++)
        {
            x(i, filter.height() - j - 1) = y(i, j);
        }
    }
    
    return x;
}

CImg<double> filter(CImg<double> input, CImg<double> filter)
{
    
    //Creates a new image with same size as the input initialized to all 0s (black)
    CImg<double> output(input, "xyzc", 0);
    
    CImg<double> filternew = flipFilter(filter);
    
    int filterWidth = filternew.width();
    int filterHeight = filternew.height();
    
    for (int i = ((filterWidth + 1) / 2) - 1; i <= input.width() - (filterWidth + 1) / 2; i++)
    {
        for (int j = ((filterHeight + 1) / 2) - 1; j <= input.height() - (filterHeight + 1) / 2; j++)
        {
            //convolution here!
            double sum = 0;
            for (int u = 0; u < filterWidth; u++)
            {
                for (int v = 0; v < filterHeight; v++)
                {
                    sum = sum + filternew(u, v) * input(i - ((filterWidth + 1) / 2) + 1 + u, j - ((filterHeight + 1) / 2) + 1 + v);
                }
            }
            
            if (sum > 255)
                sum = 255;
            
            if (sum < 0)
                sum = 0;
            
            output(i, j) = sum;
        }
    }
    
    return fixBorder(output,filterWidth);
}

CImg<double> meanFilter(CImg<double> input, int filterSize)
{
    
    if (filterSize > input.width() || filterSize > input.height())
    {
        cout << "ERROR: filter size cannot be greater than image!" << endl;
        exit(1);
    }
    
    //Creates a new grayscale image (just a matrix) to be our filter
    CImg<double> H(filterSize, filterSize, 1, 1);
    
    
    for (int i = 0; i < filterSize; i++)
    {
        for (int j = 0; j < filterSize; j++)
        {
            H(i, j) = ((double) 1) / (filterSize * filterSize);
            //H(i, j) = 1;
        }
    }
    
    //Convole with filter and return
    return filter(input, H);
    
}

CImg<double> gaussianFilter(CImg<double> input, double sigma)
{
    
    //FIXME Determine filter size, see Part 4 2b for how
    int filterSize = 3 * sigma;
    double sigmasq = sigma * sigma;
    int af = ((filterSize + 1) / 2) - 1; //adjustementFactor
    double twoPiSigmaSq = 2.0 * M_PI * sigmasq;
 
    //Creates a new grayscale image (just a matrix) to be our filter 
    CImg<double> H(filterSize, filterSize, 1, 1);
    

    for (int i = 0; i < filterSize; i++)
    {
        for (int j = 0; j < filterSize; j++)
        {
            int newi = i - af;
            int newj = j - af;
            H(i, j) = (1.0 / twoPiSigmaSq) * exp((newi * newi + newj * newj) / (-2 * sigmasq));
		}
    }
    
    H.save((char*) "Gaussiankernel.png");
    
    //Convole with filter and return
//return filter(input, H);
  return input.blur(2.5);  
}


void sort(double a[], int p, int r) //very time consuming bubble sort
{
    for (int i = 0; i <= r; i++)
    {
        for (int j = i; j <= r; j++)
        {
            if (a[i] > a[j])
            {
                int temp = a[i];
                a[i] = a[j];
                a[j] = temp;
            }
        }
    }
}

int getMedian(CImg<double> kernel)
{
    int n = kernel.width() * kernel.height();
    double a[n];  
    int temp = 0;
    for (int i = 0; i < kernel.width(); i++)
    {
        for (int j = 0; j < kernel.height(); j++)
        {
            a[temp] = kernel(i, j);
            temp++;
        }
    }
    
    //sort array now!
    sort(a, 0, n - 1);
    if (n % 2 == 0)
    {
        double x = floor((a[n / 2] + a[(n / 2) - 1]) / 2);
        return x;
    }
    else
    {
        double x = floor(a[(n - 1) / 2]);
        return x;
    }
}

CImg<double> medianFilter(CImg<double> input, int size)
{
    
    //Creates a new image with same size as the input initialized to all 0s (black)
    CImg<double> output(input, "xyzc", 0);
    CImg<double> H(size, size, 1, 1);
   
    int filterWidth = size;
    int filterHeight = size;
    
    for (int i = ((filterWidth + 1) / 2) - 1; i <= input.width() - (filterWidth + 1) / 2; i++)
    {
        for (int j = ((filterHeight + 1) / 2) - 1; j <= input.height() - (filterHeight + 1) / 2; j++)
        {
            double sum = 0;
            for (int u = 0; u < filterWidth; u++)
            {
                for (int v = 0; v < filterHeight; v++)
                {
                    H(u, v) = input(i - ((filterWidth + 1) / 2) + 1 + u, j - ((filterHeight + 1) / 2) + 1 + v);
                }
            }
            
            output(i, j) = getMedian(H);
        }
    }
    
    return fixBorder(output,filterWidth);
}

//Part 5: Quantitative Analysis

double meanSquareError(CImg<double> input, CImg<double> output)
{
	int size = input.width();
	double error = 0;
	for (int i = 0; i < input.width(); i++)
	{
		for (int j = 0; j < input.height(); j++)
		{
			error += (input(i,j)-output(i,j))*(input(i,j)-output(i,j));
		}
	}

	error = error/size;
	return error;
}
//Part 6 : Speedup

//6.1a Mean Filter
CImg<double> seperableFilter(CImg<double> input, CImg<double> filter)
{
    
    //Creates a new image with same size as the input initialized to all 0s (black)
    CImg<double> output(input, "xyzc", 0);
    
    int filterWidth = filter.width();
    int filterHeight = filter.height();
    
    for (int i = ((filterWidth + 1) / 2) - 1; i <= input.width() - (filterWidth + 1) / 2; i++)
    {
        for (int j = ((filterHeight + 1) / 2) - 1; j <= input.height() - (filterHeight + 1) / 2; j++)
        {
            //convolution here!
            double sum = 0;
            for (int u = 0; u < filterWidth; u++)
            {
                for (int v = 0; v < filterHeight; v++)
                {
                    sum = sum + filter(u, v) * input(i - ((filterWidth + 1) / 2) + 1 + u, j - ((filterHeight + 1) / 2) + 1 + v);
                }
            }
            
            if (sum > 255)
                sum = 255;
            
            if (sum < 0)
                sum = 0;
            
            output(i, j) = sum;
        }
    }
    
    return output;
}

//6.1 (Gaussian)
CImg<double> separableKCGaussian(CImg<double> input)
{
    double sigma = 3.0;
    double sigmaRoot2pi = sigma * sqrt(2.0 * M_PI);
    double minusTwoSigmaSq = -2.0 * sigma * sigma;
    
    int filterSize = 9;
    
    CImg<double> filter1(filterSize, 1, 1, 1);
    CImg<double> filter2(1, filterSize, 1, 1);
    CImg<double> firstPassImage(input, "xyzc", 0);
    CImg<double> output(input, "xyzc", 0);
    
    for (int i = 0; i < filterSize; i++)
    {
        double newi = i - (((filterSize + 1) / 2) - 1);
        filter1(i, 0) = (1.0 / sigmaRoot2pi) * exp((newi * newi) / (minusTwoSigmaSq));
    }
    
    for (int j = 0; j < filterSize; j++)
    {
        double newj = j - (((filterSize + 1) / 2) - 1);
        filter2(0, j) = (1.0 / sigmaRoot2pi) * exp((newj * newj) / (minusTwoSigmaSq));
    }
    
    return fixBorder(seperableFilter(seperableFilter(input, filter1), filter2),filterSize);
}

CImg<double> separableKernelConvolution(CImg<double> input, int filterSize)
{
	//Creates a new grayscale image (just a matrix) to be our filter
	CImg<double> filter1(filterSize, 1, 1, 1);
	CImg<double> filter2(1, filterSize, 1, 1);
	CImg<double> firstPassImage(input, "xyzc", 0);
	CImg<double> output(input, "xyzc", 0);
	
	if (filterSize > input.width() || filterSize > input.height())
	{
		cout << "ERROR: filter size cannot be greater than image!" << endl;
		exit(1);
	}
	
	//filter1 values

	for (int i = 0; i < filterSize; i++)
	{
		filter1(i, 0) = ((double) 1) / (filterSize);
	}

	//Convole with filter and return
	firstPassImage = seperableFilter(input, filter1);
	
	//filter2 values

	for (int i = 0; i < filterSize; i++)
	{
		filter2(0, i) = ((double) 1) / (filterSize);
	}
		
	//Convole with filter and return
	return fixBorder(seperableFilter(firstPassImage, filter2),filterSize);
	
}

//6.2
CImg<double> dynamicFilter1(CImg<double> input, int filterSize)
{
    CImg<double> output(input, "xyzc", 0);
    
    int filterWidth = filterSize;
    int filterHeight = 1;
    double kernelSum = 0.0;
    double sum;
    
    for (int k = 0; k < filterWidth; k++)
    {
        kernelSum = kernelSum + (input(k, 0) / filterSize);
    }
    
    for (int j = ((filterHeight + 1) / 2) - 1; j <= input.height() - ((filterHeight + 1) / 2); j++)
    {
        for (int i = ((filterWidth + 1) / 2) - 1; i <= input.width() - ((filterWidth + 1) / 2); i++)
        {
            sum = 0;
            sum = kernelSum + (input(i + (filterWidth + 1) / 2, j) / filterSize) - (input(i - ((filterWidth + 1) / 2) + 1, j) / filterSize);
            
            if (sum > 255)
                sum = 255;
            
            if (sum < 0)
                sum = 0;
            
            output(i + 1, j) = sum;
            kernelSum = sum;
        }
        
        kernelSum = 0.0;
        for (int k = 0; k < filterWidth; k++)
        {
            kernelSum = kernelSum + (input(k, j) / filterSize);
        }
    }
    return output;
}

CImg<double> dynamicFilter2(CImg<double> input, int filterSize)
{
    CImg<double> output(input, "xyzc", 0);
    
    int filterWidth = 1;
    int filterHeight = filterSize;
    double kernelSum = 0.0;
    double sum;
    
    for (int k = 0; k < filterHeight; k++)
    {
        kernelSum = kernelSum + (input(0, k) / filterSize);
    }
    
    for (int i = ((filterWidth + 1) / 2) - 1; i <= input.width() - ((filterWidth + 1) / 2); i++)
    {
        for (int j = ((filterHeight + 1) / 2) - 1; j <= input.height() - ((filterHeight + 1) / 2); j++)
        {
            sum = 0;
            sum = kernelSum + (input(i, j + ((filterHeight + 1) / 2)) / filterSize) - (input(i, j - ((filterHeight + 1) / 2) + 1) / filterSize);
            
            if (sum > 255)
                sum = 255;
            
            if (sum < 0)
                sum = 0;
            
            output(i, j + 1) = sum;
            kernelSum = sum;
        }
        
        kernelSum = 0.0;
        for (int k = 0; k < filterHeight; k++)
        {
            kernelSum = kernelSum + (input(i, k) / filterSize);
        }
    }
    
    return output;
}

CImg<double> dynamicProgrammingConvolution(CImg<double> input, int filterSize)
{
    
    CImg<double> firstPassImage(input, "xyzc", 0);
    CImg<double> output(input, "xyzc", 0);
    
    if (filterSize > input.width() || filterSize > input.height())
    {
        cout << "ERROR: filter size cannot be greater than image!" << endl;
        exit(1);
    }
    if(filterSize==1)
	cout<<"\n please enter any value greater than 1";
	
    //Convole with filter and return
    firstPassImage = dynamicFilter1(input, filterSize);
    //Convole with filter and return
    output = dynamicFilter2(firstPassImage, filterSize);
	return fixBorder(output,filterSize*2);
}

//6.3 
CImg<double> cascadedGuassianFilter(CImg<double> input, int filterSize)
{
	CImg<double> output(input, "xyzc", 0);
	CImg<double> firstPassImage(input, "xyzc", 0);
	CImg<double> image(input, "xyzc", 0);
    
	if (filterSize > input.width() || filterSize > input.height())
    {
        cout << "ERROR: filter size cannot be greater than image!" << endl;
        exit(1);
    }
    
    //Convole with filter and return
    firstPassImage = dynamicProgrammingConvolution(input, filterSize);
    
	for(int i=0;i<3;i++)
	{
		image = dynamicProgrammingConvolution(firstPassImage, filterSize);
		firstPassImage = image;
	}
	
	return fixBorder(firstPassImage,filterSize);
}