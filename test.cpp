#include <iostream>
#include <stack>
using namespace std;
class Matrix
{
public:
    int nRow, nCol;
    float matrix[4][4];

    Matrix(int rows, int cols)
    {
        nRow = rows;
        nCol = cols;
    }
    static Matrix identityMatrix(int n)
    {
        Matrix m(n,n);
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (i == j)
                    m.matrix[i][j] = 1;
                else
                    m.matrix[i][j] = 0;
            }
        }
        return m;
    }

    Matrix operator+ (Matrix mat)
    {
        Matrix m(nRow, nCol);
        for (int i = 0; i < nRow; i++)
        {
            for (int j = 0; j < nCol; j++)
            {
                m.matrix[i][j] = matrix[i][j] + mat.matrix[i][j];
            }
        }
        return m;
    }

    Matrix operator- (Matrix mat)
    {
        Matrix m(nRow, nCol);
        for (int i = 0; i < nRow; i++)
        {
            for (int j = 0; j < nCol; j++)
            {
                m.matrix[i][j] = matrix[i][j] - mat.matrix[i][j];
            }
        }
        return m;
    }

    Matrix operator* (double d)
    {
        Matrix m(nRow, nCol);
        for (int i = 0; i < nRow; i++)
        {
            for (int j = 0; j < nCol; j++)
            {
                m.matrix[i][j] = d * matrix[i][j];
            }
        }
        return m;
    }

    Matrix operator* (Matrix mat)
    {
        Matrix m(nRow, mat.nCol);
        for (int i = 0; i < m.nRow; i++)
        {
            for (int j = 0; j < m.nCol; j++)
            {
                double temp = 0;
                for (int k = 0; k < nCol; k++)
                {
                    temp += matrix[i][k] * mat.matrix[k][j];
                }
                m.matrix[i][j] = temp;
            }
        }
        return m;
    }

    void print()
    {
        for (int i = 0; i < nRow; i++)
        {
            for (int j = 0; j < nCol; j++)
            {
                cout << matrix[i][j] << "   ";
            }
            cout << endl;
        }
    }
};

int main(void)
{
    stack<string> s;
    stack<Matrix> sm;
    s.push("Hello");
    s.push("Hii");
    string str = s.top();
    cout << str;
    cout << s.size();
    cout << "Hello World." << endl;
    cout << "Bye" << endl;

    Matrix m = Matrix::identityMatrix(4);
    m.print();
    Matrix m1 = Matrix::identityMatrix(4);
    (m*(m+m1)).print();

    sm.push(m);
    sm.push(m1);
    Matrix a = sm.top();
    (a+sm.top()).print();
    cout<< sm.size();
    return 0;
}

