/*
 *  Copyright (c), 2017, Adrien Devresse
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 *
 */
#include <iostream>
#include <string>
#include <vector>

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>

const std::string FILE_NAME("create_extensible_dataset_example.h5");
const std::string DATASET_NAME("dset");

// Create a dataset name "dset" of double 4x6
//
int main(void) {
    using namespace HighFive;
    try {
        // Create a new file using the default property lists.
        File file(FILE_NAME, File::ReadWrite | File::Create | File::Truncate);

        // Create a dataspace with initial shape and max shape
        DataSpace dataspace = DataSpace({4, 5}, {17, DataSpace::UNLIMITED});

        // Use chunking
        DataSetCreateProps props;
        props.add(Chunking(std::vector<hsize_t>{2, 2}));

        // Create the dataset
        DataSet dataset = file.createDataSet(DATASET_NAME, dataspace,
                                             AtomicType<double>(), props);

        // Write into the initial part of the dataset
        double t1[3][1] = {{2.0}, {2.0}, {4.0}};
        dataset.select({0, 0}, {3, 1}).write(t1);

        // Resize the dataset to a larger size
        dataset.resize({4, 6});

        // Write into the new part of the dataset
        double t2[1][3] = {{4.0, 8.0, 6.0}};
        dataset.select({3, 3}, {1, 3}).write(t2);

        std::cout<<"Now my turn!"<<std::endl;
        dataset.resize({6, 8});
        double t3[2][3] = {{41.0, 81.0, 61.0},{41.0, 81.0, 61.0}};
        dataset.select({4, 2}, {2, 3}).write(t3);

        // now we read it back
        std::vector<std::vector<double>> result;
        dataset.read(result);


        for (auto row : result) {
            for (auto col : row)
                std::cout << " " << col;
            std::cout << std::endl;
        }

        // Now let's try adding a second dataset in the same file

        // Create the dataset
        DataSet dataset2 = file.createDataSet("Whatever", dataspace,
                                             AtomicType<double>(), props);


        double t4[2][3]={{0, 0, 1.0},{1.0, 1.0, 1.0}};
        // Write into the initial part of the dataset
        dataset2.select({0, 0}, {3, 1}).write(t1);

        // Resize the dataset to a larger size
        dataset2.resize({4, 6});

        // Write into the new part of the dataset
        dataset2.select({3, 3}, {1, 3}).write(t2);

        std::cout<<"Now my turn!"<<std::endl;
        dataset2.resize({6, 8});
        dataset2.select({4, 2}, {2, 3}).write(t4);

        // now we read it back
        dataset2.read(result);
        for (auto row : result) {
            for (auto col : row)
            std::cout << " " << col;
            std::cout << std::endl;
            }
            std::cout << "and now we are done "<< std::endl;


    } catch (const Exception& err) {
        // catch and print any HDF5 error
        std::cerr << err.what() << std::endl;
    }

    return 0; // successfully terminated
}
