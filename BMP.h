#pragma once
#include <stdio.h>
#include <vector>
using namespace std;

/**
 * @brief BMP类（存放BMP相关信息）
 * 
 */
class BMP {
public:
    // 图像类型
    unsigned short bfType;

    // 图像头结构体
    struct BitmapFileHeader {
        // unsigned short bfType;        // 19778，必须是BM字符串，对应的十六进制为0x4d42,十进制为19778，否则不是bmp格式文件
        unsigned int   bfSize;        // 文件大小 以字节为单位(2-5字节)
        unsigned short bfReserved1;   // 保留，必须设置为0 (6-7字节)
        unsigned short bfReserved2;   // 保留，必须设置为0 (8-9字节)
        unsigned int   bfOffBits;     // 从文件头到像素数据的偏移  (10-13字节)
    } bitmapFileHeader;

    // 图像信息头结构体
    struct BitmapInfoHeader {
        unsigned int    biSize;          // 此结构体的大小 (14-17字节)
        unsigned int    biWidth;         // 图像的宽  (18-21字节)
        unsigned int    biHeight;        // 图像的高  (22-25字节)
        unsigned short  biPlanes;        // 表示bmp图片的平面属，显然显示器只有一个平面，所以恒等于1 (26-27字节)
        unsigned short  biBitCount;      // 一像素所占的位数，一般为24   (28-29字节)
        unsigned int    biCompression;   // 说明图象数据压缩的类型，0为不压缩。 (30-33字节)
        unsigned int    biSizeImage;     // 像素数据所占大小, 这个值应该等于上面文件头结构中bfSize-bfOffBits (34-37字节)
        unsigned int    biXPelsPerMeter; // 说明水平分辨率，用象素/米表示。一般为0 (38-41字节)
        unsigned int    biYPelsPerMeter; // 说明垂直分辨率，用象素/米表示。一般为0 (42-45字节)
        unsigned int    biClrUsed;       // 说明位图实际使用的彩色表中的颜色索引数（设为0的话，则说明使用所有调色板项）。 (46-49字节)
        unsigned int    biClrImportant;  // 说明对图象显示有重要影响的颜色索引的数目，如果是0，表示都重要。(50-53字节)
    } bitmapInfoHeader;

    // 24位图像素信息结构体(调色板)
    struct ColorInfo {
        unsigned char rgbBlue;   //该颜色的蓝色分量  (值范围为0-255)
        unsigned char rgbGreen;  //该颜色的绿色分量  (值范围为0-255)
        unsigned char rgbRed;    //该颜色的红色分量  (值范围为0-255)
        unsigned char rgbReserved;// 保留，必须为0
    };
    vector<ColorInfo> colorInfo;

    // 图像像素信息
    vector<unsigned char> pixelInfo;

    /**
     * @brief 通过文件路径载入BMP
     * 
     * @param filePath 文件路径
     */
    BMP(const char* filePath) {
        FILE* fp;
        fopen_s(&fp, filePath, "rb");
        if (fp == NULL) {
            printf("打开图片失败!!\n");
            return;
        }
        // 读取图像类型
        fread(&bfType, 2, 1, fp);
        // 读取图像头结构体
        fread(&bitmapFileHeader, sizeof(bitmapFileHeader), 1, fp);
        // 读取图像信息头结构
        fread(&bitmapInfoHeader, sizeof(bitmapInfoHeader), 1, fp);
        // 读取调色板
        colorInfo.resize(1 << bitmapInfoHeader.biBitCount);
        fread(&colorInfo[0], sizeof(ColorInfo), colorInfo.size(), fp);
        // 读取图像像素信息
        pixelInfo.resize(bitmapInfoHeader.biWidth * bitmapInfoHeader.biHeight);
        fread(&pixelInfo[0], 1, pixelInfo.size(), fp);
        fclose(fp);
    }

    /**
     * @brief 打印BMP文件头
     * 
     */
    void showBitmapFileHeader() {
        printf("BitmapFileHeader{\n");
        printf("\tbfType:\t\t%d\n", bfType);
        printf("\tbfSize:\t\t%d\n", bitmapFileHeader.bfSize);
        printf("\tbfReserved1:\t%d\n", bitmapFileHeader.bfReserved1);
        printf("\tbfReserved2:\t%d\n", bitmapFileHeader.bfReserved2);
        printf("\tbfOffBits:\t%d\n", bitmapFileHeader.bfOffBits);
        printf("}\n");
    }

    /**
     * @brief 打印BMP信息头
     * 
     */
    void showBitmapInfoHeader() {
        printf("BitmapInfoHeader{\n");
        printf("\tbiSize:\t\t%d\n", bitmapInfoHeader.biSize);
        printf("\tbiWidth:\t%d\n", bitmapInfoHeader.biWidth);
        printf("\tbiHeight:\t%d\n", bitmapInfoHeader.biHeight);
        printf("\tbiPlanes:\t%d\n", bitmapInfoHeader.biPlanes);
        printf("\tbiBitCount:\t%d\n", bitmapInfoHeader.biBitCount);
        printf("\tbiCompression:\t%d\n", bitmapInfoHeader.biCompression);
        printf("\tbiSizeImage:\t%d\n", bitmapInfoHeader.biSizeImage);
        printf("\tbiXPelsPerMeter:%d\n", bitmapInfoHeader.biXPelsPerMeter);
        printf("\tbiYPelsPerMeter:%d\n", bitmapInfoHeader.biYPelsPerMeter);
        printf("\tbiClrUsed:\t%d\n", bitmapInfoHeader.biClrUsed);
        printf("\tbiClrImportant:\t%d\n", bitmapInfoHeader.biClrImportant);
        printf("}\n");
    }

    /**
     * @brief 打印BMP颜色盘
     * 
     */
    void showColorInfo() {
        printf("ColorInfo{\n");
        for (int i = 0; i < (1 << bitmapInfoHeader.biBitCount); i++) {
            printf("\t%0.2X %0.2X %0.2X %0.2X\n", colorInfo[i].rgbBlue, colorInfo[i].rgbGreen, colorInfo[i].rgbRed, colorInfo[i].rgbReserved);
        }
        printf("}\n");
    }

    /**
     * @brief 打印BMP像素信息
     * 
     */
    void showPixelInfo() {
        printf("PixelInfo{\n");
        for (int i = 0; i < bitmapInfoHeader.biWidth * bitmapInfoHeader.biHeight; i++) {
            printf("%0.2X ", pixelInfo[i]);
        }
        printf("\n}\n");
    }
};