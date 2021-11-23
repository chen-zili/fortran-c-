#include <iostream>

using namespace std;

#ifdef __cplusplus
extern "C" {
#endif

int temp_(int * nptr);

#ifdef __cplusplus
}
#endif

int temp_(int * nptr)
{
    cout << "!" << endl;
    cout << *nptr << endl;

    return 0;
}