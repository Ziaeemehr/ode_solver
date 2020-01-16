#include <iostream> // cout, endl
#include <vector>
#include <string>
#include <assert.h>

// using namespace std;
using std::cout;
using std::endl;
using std::string;
using std::vector;

class Myclass
{
private:
    const int a;

public:
    Myclass(const int a0, int b0) : a(a)
    {
        b = b0;
    }

    int b = 2;

    void setVariable()
    {
        b = 5;
    }
};

int main(int argc, char *argv[])
{
    Myclass myclass(1, 2);
    cout << myclass.b << endl;

    return 0;
}
