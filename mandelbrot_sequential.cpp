#include <iostream>
#include <vector>
#include <complex>
#include <fstream>
#include <chrono>
#include <cstdlib>

using namespace std;

int mandelbrot(complex<double> c, double power, int max_iter) {
    complex<double> z(0, 0);
    for (int i = 0; i < max_iter; ++i) {
        z = pow(z, power) + c;
        if (abs(z) > 2.0) {
            return i;
        }
    }
    return 0;
}

void write_output(const vector<int> &image, double elapsed_time, int width, int height,
                  double xmin, double ymin, double xmax, double ymax,
                  int iterations, double power) {
    ofstream out("mandelbrot_seq.ppm");
    out << "P3\n" << width << " " << height << "\n255\n";

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            int value = image[y * width + x];
            int r = (value * 5) % 256;
            int g = (value * 7) % 256;
            int b = (value * 11) % 256;
            out << r << " " << g << " " << b << " ";
        }
        out << "\n";
    }
    out.close();

    ofstream meta("mandelbrot_seq_metadata.txt");
    meta << "Dimensions: " << width << "x" << height << "\n";
    meta << "Region: x=[" << xmin << "," << xmax << "], y=[" << ymin << "," << ymax << "]\n";
    meta << "Iterations: " << iterations << "\n";
    meta << "Power: " << power << "\n";
    meta << "Time: " << elapsed_time << " seconds\n";
    meta.close();
}

int main(int argc, char **argv) {
    int width = 1000;
    int height = 1000;
    if (argc >= 2) {
        width = atoi(argv[1]);
        height = atoi(argv[2]);
    }
    
    double xmin = -2.0;
    double ymin = -1.5;
    double xmax = 1.0;
    double ymax = 1.5;
    int iterations = 100;
    double power = 2.0;

    cout << "Job Started" << endl;
    cout << "Sequential algorithm" << endl;
    cout << "Calculating Mandelbrot set (" << width << "x" << height << ")" << endl;

    vector<int> image(width * height);

    auto start_time = chrono::high_resolution_clock::now();

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            double re = xmin + (double(x) / width) * (xmax - xmin);
            double im = ymin + (double(y) / height) * (ymax - ymin);
            complex<double> c(re, im);

            image[y * width + x] = mandelbrot(c, power, iterations);
        }
    }

    auto end_time = chrono::high_resolution_clock::now();
    double elapsed_time = chrono::duration<double>(end_time - start_time).count();

    cout << "Finished in " << elapsed_time << " seconds." << endl;
    write_output(image, elapsed_time, width, height,
                 xmin, ymin, xmax, ymax, iterations, power);

    return 0;
}
