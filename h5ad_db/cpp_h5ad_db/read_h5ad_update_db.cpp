#include <H5Cpp.h>
#include <iostream>
#include <string>

void attr_op(H5::H5Location &loc, const std::string attr_name,
             void *operator_data) {
  std::cout << attr_name << std::endl;
}

int main() {
    filename = "../h5ads/ts_kidney.h5ad";
    dataset_name = "obs";
  // these are defined somewhere
  std::string file_name, dataset_name;


  H5::H5File file{file_name, H5F_ACC_RDONLY};
  auto dataset = file.openDataSet(dataset_name);

  dataset.iterateAttrs(attr_op);
}
