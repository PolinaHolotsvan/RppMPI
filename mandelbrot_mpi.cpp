#include <iostream>
#include <vector>
#include <complex>
#include <fstream>
#include <mpi.h>

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

void write_output(const vector<int>& image, double elapsed_time, int width, int height,
                  double xmin, double ymin, double xmax, double ymax,
                  int iterations, double power, int num_processes) {
    ofstream out("mandelbrot_mpi.ppm");
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

    ofstream meta("mandelbrot_mpi_metadata.txt");
    meta << "Dimensions: " << width << "x" << height << "\n";
    meta << "Region: x=[" << xmin << "," << xmax << "], y=[" << ymin << "," << ymax << "]\n";
    meta << "Iterations: " << iterations << "\n";
    meta << "Power: " << power << "\n";
    meta << "Time: " << elapsed_time << " seconds\n";
    meta << "Processes used: " << num_processes << "\n";
    meta.close();
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int width = 1000;
    int height = 1000;
    if (argc > 2) {
        width = stoi(argv[1]);
        height = stoi(argv[2]);
    }

    const double xmin = -2.0;
    const double ymin = -1.5;
    const double xmax = 1.0;
    const double ymax = 1.5;
    const int iterations = 100;
    const double power = 2.0;

    vector<int> global_image(width * height);
    double start_time = MPI_Wtime();

    if (rank == 0) {
        for (int p = 1; p < size; ++p) {
            int rows_per_process = height / (size - 1);
            int start_row = (p - 1) * rows_per_process;
            int end_row = (p == size - 1) ? height : start_row + rows_per_process;
            int num_rows = end_row - start_row;

            MPI_Send(&start_row, 1, MPI_INT, p, 0, MPI_COMM_WORLD);
            MPI_Send(&end_row, 1, MPI_INT, p, 0, MPI_COMM_WORLD);
        }

        for (int p = 1; p < size; ++p) {
            int rows_per_process = height / (size - 1);
            int start_row = (p - 1) * rows_per_process;
            int end_row = (p == size - 1) ? height : start_row + rows_per_process;
            int num_rows = end_row - start_row;

            vector<int> buffer(width * num_rows);
            MPI_Recv(buffer.data(), width * num_rows, MPI_INT, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            copy(buffer.begin(), buffer.end(), global_image.begin() + start_row * width);
        }

        double end_time = MPI_Wtime();
        double elapsed_time = end_time - start_time;

        cout << "Calculated Mandelbrot set (" << width << "x" << height << ")" << endl;
        cout << "Calculation completed in " << elapsed_time << " seconds." << endl;
        cout << "Using " << size << " processes" << endl;

        write_output(global_image, elapsed_time, width, height, xmin, ymin, xmax, ymax, iterations, power, size - 1);
    } else {
        int start_row, end_row;
        MPI_Recv(&start_row, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&end_row, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        int local_height = end_row - start_row;

        vector<int> local_image(width * local_height);
        for (int y = 0; y < local_height; ++y) {
            int global_y = start_row + y;
            for (int x = 0; x < width; ++x) {
                double re = xmin + (double(x) / width) * (xmax - xmin);
                double im = ymin + (double(global_y) / height) * (ymax - ymin);
                complex<double> c(re, im);
                local_image[y * width + x] = mandelbrot(c, power, iterations);
            }
        }

        MPI_Send(local_image.data(), width * local_height, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}
