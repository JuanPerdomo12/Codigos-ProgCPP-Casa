    #include <iostream>
    #include <string>
    #include <vector>

    void fill_randomly(std::vector <double> & vec){
        for (size_t i = 0; i < vec.size(); ++i) {
            vec[i] = static_cast<double>(rand()) / RAND_MAX * 100.0; // Random values between 0 and 100
        }
    }

    void print(const std::vector<double>& vec) {
    for (double x : vec) {
        std::cout << x << " ";
    }
    std::cout << "\n";
    }
    
    int main(int argc, char **argv) {
      // read variables
      if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " N" << std::endl;
        std::cerr << "N: size of the vector" << std::endl;
        return 1;
      }
      const int N = std::stoi(argv[1]);
    
      // set the vector
      std::vector<double> data;
      data.resize(N);
    
      // fill the vector
      fill_randomly(data);
    
      // print the vector
      print(data);
    
      return 0;
    }