#ifndef DNest4_TestCase_Data
#define DNest4_TestCase_Data

#include <vector>

class Data
{
  private:
    // image dimensions
    int ni, nj; // loaded from metadata_file

    // image boundaries
    double x_min, x_max;
    double y_min, y_max;
    double r_max;

    // pixel widths
    double dx, dy;

    // pixel center coordinates
    std::vector<std::vector<double>> x_rays;
    std::vector<std::vector<double>> y_rays;

    // image data
    std::vector<std::vector<double>> image; // loaded from image_file

    // image mask
    std::vector<std::vector<double>> mask;

    static Data instance;

  public:
    Data();

    void load(const char* metadata_file, const char* image_file);

    int get_ni() const { return ni; }
    int get_nj() const { return nj; }

    double get_x_min() const { return x_min; }
    double get_x_max() const { return x_max; }
    double get_y_min() const { return y_min; }
    double get_y_max() const { return y_max; }
    double get_r_max() const { return r_max; }

    double get_dx() const { return dx; }
    double get_dy() const { return dy; }
    
    const std::vector<std::vector<double>>& get_x_rays() const { return x_rays; }
    const std::vector<std::vector<double>>& get_y_rays() const { return y_rays; }
    
    const std::vector<std::vector<double>>& get_image() const { return image; }
    const std::vector<std::vector<double>>& get_mask() const { return mask; }

    static Data& get_instance() { return instance; }
};

#endif
