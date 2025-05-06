#include <iostream>
#include <limits>

int main(void) {
    std::cout << "Double min " << std::numeric_limits<double>::min() << '\n';
    std::cout << "Double max " << std::numeric_limits<double>::max() << '\n';
    std::cout << "Float min " << std::numeric_limits<float>::min() << '\n';
    std::cout << "Float max " << std::numeric_limits<float>::max() << '\n';
    std::cout << "Double lowest " << std::numeric_limits<double>::lowest() << '\n';
    std::cout << "Float lowest " << std::numeric_limits<float>::lowest() << '\n';
    std::cout << "Int min " << std::numeric_limits<int>::min() << '\n';
    std::cout << "Int max " << std::numeric_limits<int>::max() << '\n';
    std::cout << "Int lowest " << std::numeric_limits<int>::lowest() << '\n';
}