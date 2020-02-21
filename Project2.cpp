#include <iostream>
#include <cassert>
#include <cmath>	// for sqrt, fabs
#include <cfloat>	// for DBL_MAX
#include <cstdlib>	// for rand, srand
#include <ctime>	// for rand seed
#include <fstream>
#include <cstdio>	// for EOF
#include <string>
#include <vector>  // for vector distance calcs

// #include <gtkmm/application.h>
// #include <gtkmm/window.h>
// #include <gtkmm/drawingarea.h>

struct point // struct are classes that are all public
{
    static int d;
    static int i;
    double *coords;
    double *q;
    int label;
    double sum;
// constructor - a constructor initializes data members
    point(){
        this->coords = new double[d];
        //this->q = new double[d];
        this->label = 0;
        this->sum = 0;
        coords = { 0 }; // all elements 0
        q = { 0 }; // all elements 0
    }
    ~point() {}// destructor




    void print()
    {
        for (size_t i = d; i > 0; i--)
            std::cout << coords[i] << " ";
        std::cout << coords[0];
    }


// a member function double dist(point &q) that takes a reference to a second point and
// returns the Euclidean distance of the current point to the point q
    double dist(point &q)
    {
        int i;
        //if (point.size() != q.size())
        //    throw std::domain_error("vector subtraction must have equal length vectors");
        for (size_t i = d; i >= 0; i--)
            std::cout << coords[i];
            std::cout << point::q[i];
            sum += pow(coords[i]-point::q[i],2);

        return sqrt(sum);
    }
};

int main() {

    std::string str1 = "To be or not to be, that is the question";
    std::string str2 = "only ";
    std::string str3 = str1.substr(6, 12);
    str1.insert(32, str2);
    str1.replace(str1.find("to be", 0), 5, "to jump");
    str1.erase(9, 4);
    std::cout << str1 << std::endl;
    for(int i = 0; i < str3.length(); i++)
        std::cout << str3[i] << std::endl;
    return 0;
}
/*
int point::d;

class cloud
{
private:
    int d;
    int n;
    int k;

    // maximum possible number of points
    int nmax;

    point *points;
    point *centers;

public:
    cloud(int _d, int _nmax, int _k)
    {
        d = _d;
        point::d = _d;
        n = 0;
        k = _k;

        nmax = _nmax;

        points = new point[nmax];
        centers = new point[k];
    }

    ~cloud()
    {
        delete[] centers;
        delete[] points;
    }

    void add_point(point &p, int label)
    {
        assert(n < nmax);

        for(int m = 0; m < d; m++)
        {
            points[n].coords[m] = p.coords[m];
        }

        points[n].label = label;

        n++;
    }

    int get_d()
    {
        return d;
    }

    int get_n()
    {
        return n;
    }

    int get_k()
    {
        return k;
    }

    point& get_point(int i)
    {
        return points[i];
    }

    point& get_center(int j)
    {
        return centers[j];
    }

    void set_center(point &p, int j)
    {
        for(int m = 0; m < d; m++)
            centers[j].coords[m] = p.coords[m];
    }

    double intracluster_variance()
    {
        // TODO
        return 0.0;
    }

    void set_voronoi_labels()
    {
        // TODO
    }

    void set_centroid_centers()
    {
        // TODO
    }

    void kmeans()
    {
        set_centroid_centers(); // trivial initialization

        // TODO
    }

    void init_forgy()
    {
        // TODO
    }

    void init_plusplus()
    {
        // TODO
    }

    void init_random_partition()
    {
        // TODO
    }
};


// test functions
void test_intracluster_variance()
{
    // tolerance for comparison of doubles
    const double eps = 0.0001;

    // dimension used for tests
    point::d = 1;

    // temporary container
    point p;

    std::cout << "Testing function intracluster_variance()...";

    // test case 1
    const double dist_onepoint_zerodist = 0.0;
    cloud onepoint_zerodist(1, 1, 1);
    p.coords[0] = 0.0;
    onepoint_zerodist.add_point(p, 0);
    assert(std::fabs(onepoint_zerodist.intracluster_variance() - dist_onepoint_zerodist) < eps);

    // test case 2
    const double dist_onepoint_posdist = 0.25;
    cloud onepoint_posdist(1, 1, 1);
    p.coords[0] = 0.5;
    onepoint_posdist.add_point(p, 0);
    assert(std::fabs(onepoint_posdist.intracluster_variance() - dist_onepoint_posdist) < eps);

    // test case 3
    const double dist_twopoints = 0.625;
    cloud twopoints(1, 2, 1);
    p.coords[0] = -1.0;
    twopoints.add_point(p, 0);
    p.coords[0] = 0.5;
    twopoints.add_point(p, 0);
    p.coords[0] = -0.5;
    twopoints.set_center(p, 0);
    assert(std::fabs(twopoints.intracluster_variance() - dist_twopoints) < eps);

    // test case 4
    const double dist_twoclusters = 6.8125;
    cloud twoclusters(1, 4, 2);
    p.coords[0] = -1.0;
    twoclusters.add_point(p, 0);
    p.coords[0] = 0.5;
    twoclusters.add_point(p, 0);
    p.coords[0] = -0.5;
    twoclusters.set_center(p, 0);
    p.coords[0] = -2.0;
    twoclusters.add_point(p, 1);
    p.coords[0] = 2.0;
    twoclusters.add_point(p, 1);
    p.coords[0] = -3.0;
    twoclusters.set_center(p, 1);
    assert(std::fabs(twoclusters.intracluster_variance() - dist_twoclusters) < eps);

    std::cout << "\t\t[OK]" << std::endl;
}

void test_kmeans()
{
    // TODO
}

void test_init_forgy()
{
    // number of random experiments
    const int K = 10000;
    // tolerance in probability
    const double delta = 0.0625;

    // dimenstion used for tests
    point::d = 1;

    // temporary container
    point p;

    std::cout << "Testing function init_forgy()...";

    const double prob_threepoints = 0.3333;
    cloud threepoints(1, 3, 1);
    p.coords[0] = 0.0;
    threepoints.add_point(p, 0);
    p.coords[0] = 1.0;
    threepoints.add_point(p, 0);
    p.coords[0] = 2.0;
    threepoints.add_point(p, 0);
    int cnt = 0;
    for(int k = 0; k < K; k++)
    {
        threepoints.init_forgy();
        if(threepoints.get_center(0).coords[0] == 1.0)
            cnt++;
    }
    assert(std::fabs(cnt/(double)K - prob_threepoints) < delta);

    std::cout << "\t\t[OK]" << std::endl;
}

void test_init_plusplus()
{
    // number of random experiments
    const int K = 10000;
    // tolerance in probability
    const double delta = 0.0625;

    // dimenstion used for tests
    point::d = 1;

    // temporary container
    point p;

    std::cout << "Testing function init_plusplus()...";

    // test case 1
    const double prob_threepoints = 0.3333;
    cloud threepoints(1, 3, 1);
    p.coords[0] = 0.0;
    threepoints.add_point(p, 0);
    p.coords[0] = 1.0;
    threepoints.add_point(p, 0);
    p.coords[0] = 2.0;
    threepoints.add_point(p, 0);
    int cnt = 0;
    for(int k = 0; k < K; k++)
    {
        threepoints.init_plusplus();
        if(threepoints.get_center(0).coords[0] == 1.0)
            cnt++;
    }
    assert(std::fabs(cnt/(double)K - prob_threepoints) < delta);

    // test case 2
    const double prob_twoclusters = 0.125;
    cloud twoclusters(1, 4, 2);
    p.coords[0] = 0.0;
    twoclusters.add_point(p, 0);
    p.coords[0] = 0.0;
    twoclusters.add_point(p, 0);
    p.coords[0] = 1.0;
    twoclusters.add_point(p, 0);
    p.coords[0] = 2.0;
    twoclusters.add_point(p, 0);
    cnt = 0;
    for(int k = 0; k < K; k++)
    {
        twoclusters.init_plusplus();
        if(twoclusters.get_center(1).coords[0] == 1.0)
            cnt++;
    }
    assert(std::fabs(cnt/(double)K - prob_twoclusters) < delta);

    std::cout << "\t\t[OK]" << std::endl;
}

void test_init_random_partition()
{
    // number of random experiments
    const int K = 10000;
    // tolerance in probability
    const double delta = 0.0625;

    // dimenstion used for tests
    point::d = 1;

    // temporary container
    point p;

    std::cout << "Testing function init_random_parition()...";

    const double prob_threepoints = 0.3333;
    cloud threepoints(1, 3, 3);
    p.coords[0] = 0.0;
    threepoints.add_point(p, 0);
    p.coords[0] = 1.0;
    threepoints.add_point(p, 0);
    p.coords[0] = 2.0;
    threepoints.add_point(p, 0);
    int cnt = 0;
    for(int k = 0; k < K; k++)
    {
        threepoints.init_random_partition();
        if(threepoints.get_point(2).label == 1)
            cnt++;
    }
    assert(std::fabs(cnt/(double)K - prob_threepoints) < delta);

    std::cout << "\t[OK]" << std::endl;
}*/
// for graphical interface
/*
// for graphical interface
class MyArea : public Gtk::DrawingArea
{
private:
    cloud *c;

public:
    MyArea(cloud *_c)
    {
        c = _c;
    }

    virtual ~MyArea() {}

protected:
    //Override default signal handler:
    bool on_draw(const Cairo::RefPtr<Cairo::Context> &cr) override;
};


int main(int argc, char **argv)
{
    srand(time(NULL));

    point::d = 4;
    // TODO: your tests for the point class here

    test_intracluster_variance();
    test_kmeans();
    // test_init_forgy();
    // test_init_plusplus();
    // test_init_random_partition();

    const int d = 4;
    const int nmax = 150;
    const int k = 3;

    // construct point cloud
    cloud c(d, nmax, k);

    // open data file
    std::ifstream is("iris.data");
    assert(is.is_open());

    // point to read into
    point p;

    // labels to cycle through
    int label = 0;

    // while not at end of file
    while(is.peek() != EOF)
    {
        // read new points
        for(int m = 0; m < d; m++)
        {
            is >> p.coords[m];
        }

        c.add_point(p, label);

        label = (label + 1) % k;

        // read ground-truth labels
        // unused in normal operation
        std::string next_name;
        is >> next_name;

        // consume \n
        is.get();
    }

    // execute k-means algorithm
    std::cout << "Intracluster variance before k-means: " << c.intracluster_variance() << std::endl;
    c.kmeans();
    std::cout << "Intracluster variance after k-means: " << c.intracluster_variance() << std::endl;


    // launch graphical interface
    Glib::RefPtr<Gtk::Application> app = Gtk::Application::create(argc, argv, "inf442.td3");

    Gtk::Window win;
    win.set_title("TD 3");
    win.set_default_size(400, 400);

    MyArea area(&c);
    win.add(area);
    area.show();

    return app->run(win);
}


bool MyArea::on_draw(const Cairo::RefPtr<Cairo::Context> &cr)
{
    Gtk::Allocation allocation = get_allocation();
    const int width = allocation.get_width();
    const int height = allocation.get_height();

    // choose which coordinates to use for 2d image
    const int x_coord = 2;
    const int y_coord = 3;

    // find min and max on each axis
    double x_min = DBL_MAX;
    double x_max = DBL_MIN;
    double y_min = DBL_MAX;
    double y_max = DBL_MIN;
    for(int i = 0; i < c->get_n(); i++)
    {
        if(c->get_point(i).coords[x_coord] < x_min)
            x_min = c->get_point(i).coords[x_coord];

        if(c->get_point(i).coords[x_coord] > x_max)
            x_max = c->get_point(i).coords[x_coord];

        if(c->get_point(i).coords[y_coord] < y_min)
            y_min = c->get_point(i).coords[y_coord];

        if(c->get_point(i).coords[y_coord] > y_max)
            y_max = c->get_point(i).coords[y_coord];
    }

    // plot all points
    for(int i = 0; i < c->get_n(); i++)
    {
        cr->save(); // save current drawing context (opaque black)
        cr->arc((c->get_point(i).coords[x_coord]-x_min)*width/(x_max-x_min), (c->get_point(i).coords[y_coord]-y_min)*height/(y_max-y_min), 10.0, 0.0, 2.0 * M_PI); // full circle

        // choose color depending on label
        switch(c->get_point(i).label)
        {
            case 0:
                cr->set_source_rgba(1.0, 0.0, 0.0, 0.6); // red, partially translucent
                break;

            case 1:
                cr->set_source_rgba(0.0, 0.0, 0.8, 0.6); // 0.8 blue, partially translucent
                break;

            case 2:
                cr->set_source_rgba(0.0, 1.0, 0.0, 0.6); // green, partially translucent
                break;

            default:
                cr->set_source_rgba(1.0, 1.0, 0.0, 0.6); // yellow, partially translucent
                break;
        }

        cr->fill_preserve();
        cr->restore();  // restore drawing context to opaque black
        cr->stroke();
    }
    return true;
}
*/
