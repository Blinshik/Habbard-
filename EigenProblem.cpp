#include <iostream>
#include </progiFizika/Habbard++/EigenProblem/include/armadillo>

using namespace std;
using namespace arma;


int main()
{
    float U = 20;
    float Jd = 5;
    float U2 = 18;
    float J = 5;
    float delta = 0.1;
    float lambda = 0.3;


    vec eigval;
    mat eigvec;


    mat A = {
        {U+2*delta,Jd,Jd,0,0,0,lambda,0,0,0,0,lambda,0,0,0},
        {Jd,U,Jd,0,0,0,0,0,0,0,0,0,0,0,lambda},
        {Jd,Jd,U,lambda,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,lambda,U2+delta-J+0.5*lambda,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,U2 + delta + 0.5 * lambda,-J,0,lambda,0,0,0,0,0,0,0},
        {0,0,0,0,-J,U2 + delta - 0.5 * lambda,0,0,0,0,0,0,0,0,0},
        {lambda,0,0,0,0,0,U2 + delta - J - 0.5 * lambda,0,0,lambda,0,0,0,0,0},
        {0,0,0,0,lambda,0,0,U2-J,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,U2+lambda,-J,0,0,0,0,0},
        {0,0,0,0,0,0,lambda,0,-J,U2-lambda,0,lambda,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,U2-J,0,lambda,0,0},
        {lambda,0,0,0,0,0,0,0,0,lambda,0,U2 + delta - J - 0.5 * lambda,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,lambda,0,U2 + delta + 0.5 * lambda,-J,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,-J,U2 + delta - 0.5 * lambda,0},
        {0,lambda,0,0,0,0,0,0,0,0,0,0,0,0,U2 + delta - J + 0.5 * lambda}
    };

    

    eig_sym(eigval, eigvec, A);
    A.print("A: ");
    eigval.print("Avalue:");
    eigvec.print("Avector:");


    return 0;
}


