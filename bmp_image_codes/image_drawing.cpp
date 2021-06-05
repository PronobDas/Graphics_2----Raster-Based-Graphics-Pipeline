/*
*****************************************
* Thanks to Siddhartha Shankar Das (Sir *                                                                        *
*****************************************
*/

#include "bitmap_image.hpp"

using namespace std;


int main(){

    bitmap_image image(500,300);
    image.set_pixel(0,100, 255, 255, 255);
    image.set_pixel(1,101, 255, 255, 255);
    image.set_pixel(2,102, 255, 255, 255);
    image.set_pixel(3,103, 255, 255, 255);

    for(int i=0;i<200;i++){
        for(int j=0;j<100;j++){
            //image.set_pixel(i,j,255,255,0);
        }
        for(int j=150;j<200;j++){
            //image.set_pixel(i,j,0,0,255);
        }
    }

    for(int i = 0; i < 119; i++)
        for (int j = 0; j < 119; j++)
            image.set_pixel(i,j,255, 255, 255);


    for(int i=400;i<450;i++){
        for(int j=50;j<150;j++){
            //image.set_pixel(i,j,255,0,0);
        }
        for(int j=200;j<300;j++){
            //image.set_pixel(i,j,0,255,255);
        }
    }

    image.save_image("test2.bmp");;

    return 0;
}
