#include <iostream>
#include <stack>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "bitmap_image.hpp"
using namespace std;

#define PI 2*acos(0.0)

class Point
{
public:
    double x, y, z, w;

    Point(){}
    Point(double a, double b, double c, double d)
    {
        x = a;
        y = b;
        z = c;
        w = d;
    }

    double distance(Point p)
    {
        return sqrt(pow((p.x-x),2) + pow((p.y-y),2) + pow((p.z-z),2) );
    }

    Point operator+( Point p)
    {
        Point temp = Point();
        temp.x = x + p.x;
        temp.y = y + p.y;
        temp.z = z + p.z;
        temp.w = w + p.w;
        return temp;
    }
    Point operator-( Point p)
    {
        Point temp = Point();
        temp.x = x - p.x;
        temp.y = y - p.y;
        temp.z = z - p.z;
        temp.w = w - p.w;
        return temp;
    }
    void print()
    {
        cout << x << "  " << y << "  " << z  << "  "<< w << endl;
    }
};


class Vector
{
public:
    double x, y, z;

    Vector(){}

    Vector(double a, double b, double c)
    {
        x = a;
        y = b;
        z = c;
    }

    void normalize()
    {
        double r = sqrt(x*x + y*y + z*z);
        x = x/r;
        y = y/r;
        z = z/r;
    }

    Vector operator+( Vector v)
    {
        Vector temp = Vector();
        temp.x = x + v.x;
        temp.y = y + v.y;
        temp.z = z + v.z;
        return temp;
    }

    Vector operator-( Vector v)
    {
        Vector temp = Vector();
        temp.x = x - v.x;
        temp.y = y - v.y;
        temp.z = z - v.z;
        return temp;
    }

    // constant multiplication
    Vector operator*( double d)
    {
        Vector temp = Vector();
        temp.x = d*x;
        temp.y = d*y;
        temp.z = d*z;
        return temp;
    }

    static double dotMul(Vector a, Vector b)
    {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }

    static Vector crossMul(Vector a, Vector b)
    {
        Vector temp = Vector();
        temp.x = a.y * b.z - b.y * a.z;
        temp.y = b.x * a.z - a.x * b.z;
        temp.z = a.x * b.y - b.x * a.y;
        return temp;
    }

    void print()
    {
        cout << x << "  " << y << "  " << z << endl;
    }
};


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


Point transformPoint( Matrix mat, Point p)
{
    Matrix m(4, 1);
    m.matrix[0][0] = p.x;
    m.matrix[1][0] = p.y;
    m.matrix[2][0] = p.z;
    m.matrix[3][0] = p.w;

    Matrix temp = mat * m;
    Point p_new(temp.matrix[0][0], temp.matrix[1][0], temp.matrix[2][0], temp.matrix[3][0]);
    return p_new;
}


Vector rodrigues( Vector x, Vector a, double ang)
{
    double angle = ang*PI/180;
    double temp1 = Vector::dotMul(a, x);
    Vector temp2 = Vector::crossMul(a, x);

    double sTheta = sin(angle);
    double cTheta = sqrt( 1 - sTheta*sTheta );

    Vector temp = x*cTheta + (a*(1-cTheta))*temp1 + temp2*sTheta;
    return temp;
}


class Color
{
public:
    double R, G, B;

    Color(){}

    Color(double r, double g, double b)
    {
        R = r;
        G = g;
        B = b;
    }

    void print()
    {
        cout << R << " " << G << " " << B << endl;
    }
};


class Triangle
{
public:
    Point points[3];
    Color color;

    Triangle(){}

    Triangle(Point p1, Point p2, Point p3, Color c)
    {
        points[0] = p1;
        points[1] = p2;
        points[2] = p3;
        color = c;
    }
};


int main(void)
{
    stack<Matrix> sm;
    Matrix m = Matrix::identityMatrix(4);
    sm.push(m);

    Matrix currentTop = sm.top();

    double eyeX, eyeY, eyeZ;
    double lookX, lookY, lookZ;
    double upX, upY, upZ;
    double fovY, aspectRatio, near, far;

    ifstream scene;
    ofstream stage1;
    scene.open("scene.txt");
    stage1.open("stage1.txt");
    stage1 << std::fixed << std::setprecision(7);

    scene >> eyeX >> eyeY >> eyeZ;
    scene >> lookX >> lookY >> lookZ;
    scene >> upX >> upY >> upZ;
    scene >> fovY >> aspectRatio >> near >> far;

    int tCount = 0;

    // stage 1: Modeling Transformation
    string cmd;
    while(true)
    {
        scene >> cmd;
        if (cmd == "triangle")
        {
            tCount++;
            //cout << "triangle" << endl;
            for (int i = 0; i < 3; i++)
            {
                double d1, d2, d3;
                scene >> d1 >> d2 >> d3;

                Point p(d1, d2, d3, 1);
                Point temp = transformPoint( currentTop, p );
                //temp.print();

                stage1 << temp.x << "  " << temp.y << "  " << temp.z << endl;
            }
            stage1 << endl;
        }
        else if (cmd == "translate")
        {
            //cout << "translate"<< endl;
            double tx, ty, tz;
            scene >> tx >> ty >> tz;
            Matrix temp = Matrix::identityMatrix(4);
            temp.matrix[0][3] = tx;
            temp.matrix[1][3] = ty;
            temp.matrix[2][3] = tz;

            currentTop = currentTop * temp;

        }

        else if (cmd == "scale")
        {
            //cout << "scale"<< endl;
            double sx, sy, sz;
            scene >> sx >> sy >> sz;
            Matrix temp = Matrix::identityMatrix(4);
            temp.matrix[0][0] = sx;
            temp.matrix[1][1] = sy;
            temp.matrix[2][2] = sz;

            currentTop = currentTop * temp;
        }

        else if (cmd == "rotate")
        {
            //cout << "rotate"<< endl;

            double angle, ax, ay, az;
            scene >> angle >> ax >> ay >> az;

            Vector a(ax, ay, az);
            a.normalize();

            Vector i(1, 0, 0);
            Vector j(0, 1, 0);
            Vector k(0, 0, 1);

            Vector c1 = rodrigues(i, a, angle);
            Vector c2 = rodrigues(j, a, angle);
            Vector c3 = rodrigues(k, a, angle);

            Matrix temp = Matrix::identityMatrix(4);
            temp.matrix[0][0] = c1.x;
            temp.matrix[1][0] = c1.y;
            temp.matrix[2][0] = c1.z;

            temp.matrix[0][1] = c2.x;
            temp.matrix[1][1] = c2.y;
            temp.matrix[2][1] = c2.z;

            temp.matrix[0][2] = c3.x;
            temp.matrix[1][2] = c3.y;
            temp.matrix[2][2] = c3.z;

            currentTop = currentTop * temp;
        }

        else if (cmd == "push")
        {
            //cout << "push"<< endl;
            sm.push(currentTop);
            //currentTop = Matrix::identityMatrix(4);
        }

        else if (cmd == "pop")
        {
            //cout << "pop"<< endl;
            currentTop = sm.top();
            sm.pop();
        }

        else if (cmd == "end")
        {
            break;
        }
    }

    scene.close();
    stage1.close();



    // stage 2: View Transformation
    Vector look(lookX, lookY, lookZ);
    Vector eye(eyeX, eyeY, eyeZ);
    Vector up(upX, upY, upZ);

    Vector l = look - eye;
    l.normalize();

    Vector r = Vector::crossMul(l, up);
    r.normalize();

    Vector u = Vector::crossMul(r, l);

    Matrix T = Matrix::identityMatrix(4);
    T.matrix[0][3] = -eyeX;
    T.matrix[1][3] = -eyeY;
    T.matrix[2][3] = -eyeZ;

    Matrix R = Matrix::identityMatrix(4);
    R.matrix[0][0] = r.x;
    R.matrix[0][1] = r.y;
    R.matrix[0][2] = r.z;

    R.matrix[1][0] = u.x;
    R.matrix[1][1] = u.y;
    R.matrix[1][2] = u.z;

    R.matrix[2][0] = -l.x;
    R.matrix[2][1] = -l.y;
    R.matrix[2][2] = -l.z;

    Matrix V = R * T;

    ifstream stageIn;
    ofstream stage2;
    stageIn.open("stage1.txt");
    stage2.open("stage2.txt");
    stage2 << std::fixed << std::setprecision(7);

    int tCountTemp = tCount;
    while( tCountTemp )
    {
        for (int i = 0;  i < 3; i++)
        {
            double d1, d2, d3;
            stageIn >> d1 >> d2 >> d3;

            Point p(d1, d2, d3, 1);
            Point temp = transformPoint( V , p );
    //        temp.print();

            stage2 << temp.x << "  " << temp.y << "  " << temp.z << endl;
        }
        stage2 << endl;
        tCountTemp--;
    }
    stageIn.close();
    stage2.close();



    // stage 3 : Projection Transformation
    double fovX = fovY * aspectRatio;
    double _t = near * tan(PI * fovY / 360);
    double _r = near * tan(PI * fovX/ 360);

    Matrix P = Matrix::identityMatrix(4);
    P.matrix[0][0] = near/_r;
    P.matrix[1][1] = near/_t;
    P.matrix[2][2] = -(far + near) / (far - near);
    P.matrix[2][3] = -(2 * far * near) / (far - near);
    P.matrix[3][2] = -1;
    P.matrix[3][3] = 0;

    ifstream stageIn2;
    ofstream stage3;
    stageIn2.open("stage2.txt");
    stage3.open("stage3.txt");
    stage3 << std::fixed << std::setprecision(7);

    tCountTemp = tCount;
    while( tCountTemp )
    {
        for (int i = 0;  i < 3; i++)
        {
            double d1, d2, d3;
            stageIn2 >> d1 >> d2 >> d3;

            Point pPoint(d1, d2, d3, 1);
            //pPoint.print();
            Point temp = transformPoint( P , pPoint );
            //P.print();
            //temp.print();

            stage3 << temp.x / temp.w << "  " << temp.y / temp.w<< "  " << temp.z / temp.w<< endl;
        }
        stage3 << endl;
        tCountTemp--;
    }
    stageIn2.close();
    stage3.close();




    // stage 4 : Clipping & scan conversion using Z-buffer algorithm

    ifstream stageIn3;
    ifstream config;
    ofstream z_buffer;
    stageIn3.open("stage3.txt");
    config.open("config.txt");
    z_buffer.open("z_buffer.txt");
    z_buffer << std::fixed << std::setprecision(6);

    int screenWidth, screenHeight;
    double leftLimit, rightLimit, bottomLimit, topLimit, frontLimit, rearLimit;

    config >> screenWidth >> screenHeight;
    config >> leftLimit >> bottomLimit;
    rightLimit = -leftLimit;
    topLimit = -bottomLimit;
    config >> frontLimit >> rearLimit;
    config.close();

    Triangle triangles[tCount];

    //storing the file inputs into Triangle data structure
    tCountTemp = tCount;
    while (tCountTemp)
    {
        Point p[3];
        for(int i = 0; i < 3; i++)
        {
            stageIn3 >> p[i].x >> p[i].y >> p[i].z;
            p[i].w = 1;
        }

        Color c(rand()%255, rand()%255, rand()%255);
        Triangle temp(p[0], p[1], p[2], c);
        triangles[tCount-tCountTemp] = temp;

        tCountTemp--;
    }


    //initializing necessary variables
    double dx = (rightLimit - leftLimit)/screenWidth;
    double dy = (topLimit - bottomLimit)/screenHeight;

    double topY = topLimit - dy/2;
    double leftX = leftLimit + dx/2;
    double bottomY = bottomLimit + dy/2;
    double rightX = rightLimit - dx/2;


    //initializing z_buffer array
    double** zBufferArray = new double*[screenHeight];
    for (int i = 0; i < screenWidth; i++)
        zBufferArray[i] = new double[screenWidth];

    for (int i = 0; i < screenHeight; i++)
    {
        for (int j = 0; j < screenWidth; j++)
            zBufferArray[i][j] = rearLimit;
    }


    //initializing image buffer and setting background as black.
    bitmap_image image(screenWidth, screenHeight);
    for(int i = 0; i < screenWidth; i++)
    {
        for (int j = 0; j < screenHeight; j++)
            image.set_pixel(i,j,0, 0, 0);
    }



    //Start of Main Algorithm

    //for every triangle
    for (int i = 0; i < tCount; i++)
    {
        Triangle triangle = triangles[i];
        double a = triangle.points[0].y;
        double b = triangle.points[1].y;
        double c = triangle.points[2].y;

        // maxX, maxY, minX, minY refers max and min value of the triangle
        double maxY=( a>=b&&a>=c)?a:(b>=a&&b>=c)?b:c;
        double minY=(a<=b&&a<=c)?a:(b<=a&&b<=c)?b:c;

        // if max, min values go outside of the viewing area, clipped them.
        if (maxY > topY)
            maxY = topY;
        if (minY < bottomY)
            minY = bottomY;


        a = triangle.points[0].x;
        b = triangle.points[1].x;
        c = triangle.points[2].x;

        double maxX=(a>=b&&a>=c)?a:(b>=a&&b>=c)?b:c;
        double minX=(a<=b&&a<=c)?a:(b<=a&&b<=c)?b:c;

        if (maxX > rightX)
            maxX = rightX;
        if (minX < leftX)
            minX = leftX;


        //calculate top scanLine and bottom scanLine
        int topRow = round((topY - maxY)/dy);
        int bottomRow = round((topY - minY)/dy);

        //for each scanLine
        for (int row = topRow; row < bottomRow; row++)
        {
            // ys--> y co-ordinate of a specific row
            double ys = topY - row * dy;

            //for finding left starting point and right starting point of a scanLine
            double xa = 1000, xb = 1000, xc = 1000, za = 1000, zb = 1000, zc = 1000;

            //calculating 3 possible intersecting point.
            if ( (triangle.points[1].y - triangle.points[0].y) != 0 )
            {
                xa = triangle.points[0].x + (ys - triangle.points[0].y)*(triangle.points[1].x - triangle.points[0].x)/(triangle.points[1].y - triangle.points[0].y);
               // za = triangle.points[0].z + (ys - triangle.points[0].y)*(triangle.points[1].z - triangle.points[0].z)/(triangle.points[1].y - triangle.points[0].y);
            }

            if ( (triangle.points[2].y - triangle.points[0].y) != 0 )
            {
                xb = triangle.points[0].x + (ys - triangle.points[0].y)*(triangle.points[2].x - triangle.points[0].x)/(triangle.points[2].y - triangle.points[0].y);
               // zb = triangle.points[0].z + (ys - triangle.points[0].y)*(triangle.points[2].z - triangle.points[0].z)/(triangle.points[2].y - triangle.points[0].y);
            }

            if ( (triangle.points[2].y - triangle.points[1].y) != 0 )
            {
                xc = triangle.points[1].x + (ys - triangle.points[1].y)*(triangle.points[2].x - triangle.points[1].x)/(triangle.points[2].y - triangle.points[1].y);
               // zc = triangle.points[1].z + (ys - triangle.points[1].y)*(triangle.points[2].z - triangle.points[1].z)/(triangle.points[2].y - triangle.points[1].y);
            }

            // for leftX and rightX value of a row(scanLine)
            double maxa;
            double mina;

            if ( xa == 1000 )
            {
                maxa = (xb >= xc)?xb:xc;
                mina = (xb <= xc)?xb:xc;
            }
            else if ( xb == 1000 )
            {
                maxa = (xa >= xc)?xa:xc;
                mina = (xa <= xc)?xa:xc;
            }
            else if ( xc == 1000 )
            {
                maxa = (xa >= xb)?xa:xb;
                mina = (xa <= xb)?xa:xb;
            }
            else if ( (xa == 1000 && xb == 1000) || (xc == 1000 && xb == 1000) || (xa == 1000 && xc == 1000) )
            {
                cout << "continue" << endl;
                continue;
            }
            else
            {
                if ( xa > maxX || xa < minX )
                {
                    maxa = (xb >= xc)?xb:xc;
                    mina = (xb < xc)?xb:xc;
                }
                else if ( xb > maxX || xb < minX  )
                {
                    maxa = (xa >= xc)?xa:xc;
                    mina = (xa < xc)?xa:xc;
                }
                else
                {
                    maxa = (xb >= xa)?xb:xa;
                    mina = (xb < xa)?xb:xa;
                }
            }


            //calculating z values for the two intersecting points
            //za --> z value of left intersection point. zb --> right intersection point.
            if ( maxa == xa)
                zb = triangle.points[0].z + (ys - triangle.points[0].y)*(triangle.points[1].z - triangle.points[0].z)/(triangle.points[1].y - triangle.points[0].y);
            else if (maxa == xb)
                zb = triangle.points[0].z + (ys - triangle.points[0].y)*(triangle.points[2].z - triangle.points[0].z)/(triangle.points[2].y - triangle.points[0].y);
            else
                zb = triangle.points[1].z + (ys - triangle.points[1].y)*(triangle.points[2].z - triangle.points[1].z)/(triangle.points[2].y - triangle.points[1].y);


            if ( mina == xa)
                za = triangle.points[0].z + (ys - triangle.points[0].y)*(triangle.points[1].z - triangle.points[0].z)/(triangle.points[1].y - triangle.points[0].y);
            else if (mina == xb)
                za = triangle.points[0].z + (ys - triangle.points[0].y)*(triangle.points[2].z - triangle.points[0].z)/(triangle.points[2].y - triangle.points[0].y);
            else
                za = triangle.points[1].z + (ys - triangle.points[1].y)*(triangle.points[2].z - triangle.points[1].z)/(triangle.points[2].y - triangle.points[1].y);


            // if max, min values go outside of the viewing area, clipped them.
            if (maxa > rightX)
                maxa = rightX;
            if (mina < leftX)
                mina = leftX;


            //calculating left and right column number.
            int leftCol = round((mina - leftX)/dx);
            int rightCol = round((maxa - leftX)/dx);

            //for each column(in which we need to set pixel value) in a scanLine.
            for (int col = leftCol; col < rightCol; col++)
            {
                double xp = leftX + col * dx;
                //updating z value in a corresponding pixel
                double zp = za + (xp - mina) * (zb - za)/(maxa - mina);

                if ( zp > frontLimit && zp < zBufferArray[row][col] )
                {
                    zBufferArray[row][col] = zp;
                    image.set_pixel(col, row, triangle.color.R, triangle.color.G, triangle.color.B);

                    z_buffer << zp << "\t" ;
                }
            }
            z_buffer << endl;
        }
    }


    image.save_image("output.bmp");

    stageIn3.close();
    z_buffer.close();

    //free the z_buffer array
    for (int i = 0; i < screenWidth; i++)
        delete[] zBufferArray[i];

    return 0;
}
