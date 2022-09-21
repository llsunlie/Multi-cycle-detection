/**
 * @file MyMulti-cycle detection.h
 * @author 201905555514 刘世豪
 * @brief 我的多圆检测
 * @version 0.1
 * @date 2022-05-23
 * 
 * @copyright Copyright (c) 2022 liushihao
 * 
 */
#pragma once
#include <set>
#include <float.h>
#include <cmath>
#include <cstring>
#include <time.h>
#include "BMP.h"
#include "Point.h"
#include "Cycle.h"

// 每次寻找一个圆最大的回合数
#define MAX_ROUNDS 20000
// 最多寻找两个有效点的回合数
#define MAX_FIND_TWO_POINTS_ROUNDS 30
// 寻找两个有效点最小距离
#define MIN_TWO_POINTS_DISTANCE 5
// 两直线角度误差
#define ANGLE_ERROR 5
// 判断三点共圆的阈值
#define JUDGE_DIRECTION_THRESHOLD 1
// 比例系数
#define PROPORTIONALITY_FACTOR 0.75
// 距离误差
#define DISTANCE_ERROR 0.9

const double PI = acos(-1.0);

class MyMulticycleDetection {
public:
	// 图像高度
	int height;
	// 图像宽度
	int width;
	// 初始图像矩阵
	vector<vector<int>> originalMatrix;
	// 最终图像矩阵
	vector<vector<int>> finalMatrix;
	// 候选点集
	set<Point> candidateSet;
	// 结果圆集合
    vector<Cycle> resultCycleSet;

	/**
	 * @brief Construct a new My Multicycle Detection object
	 * 
	 * @param bmp 
	 */
	MyMulticycleDetection(const BMP &bmp) {
		// 初始化高宽
		this->height = bmp.bitmapInfoHeader.biHeight;
		this->width = bmp.bitmapInfoHeader.biWidth;
		// 初始化初始图矩阵和最终图矩阵
		this->originalMatrix.resize(height);
		this->finalMatrix.resize(height);
		for (int i = 0; i < height; i++) {
			this->originalMatrix[height - i - 1].resize(width);
			this->finalMatrix[i].resize(width);
			for (int j = 0; j < width; j++) {
				this->originalMatrix[height - i - 1][j] =
					((int)bmp.pixelInfo[i * width + j] == 0);
			}
		}
		// 初始化候选点集
		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				if (originalMatrix[i][j]) {
					Point tmp(i, j);
					candidateSet.insert(tmp);
				}
			}
		}
	}

	/**
	 * @brief 将结果输出到文件
	 * 
	 * @param bmp 原图像 bmp 图像
	 * @param filePath 保存的文件路径
	 */
	void finalBMPToFile(const BMP& bmp, const char* filePath) {
		FILE* fp;
		fopen_s(&fp, filePath, "wb");
		// 输出图像类型
		fwrite(&bmp.bfType, 2, 1, fp);
		// 输出图像头结构体
		fwrite(&bmp.bitmapFileHeader, sizeof(bmp.bitmapFileHeader), 1, fp);
		// 输出图像信息头结构
		fwrite(&bmp.bitmapInfoHeader, sizeof(bmp.bitmapInfoHeader), 1, fp);
		// 输出调色板
		fwrite(&bmp.colorInfo[0], sizeof(BMP::ColorInfo), bmp.colorInfo.size(), fp);
		// 输出图像像素信息
		for (int i = height - 1; i >= 0; i--) {
			for (int j = 0; j < width; j++) {
				unsigned char c = (finalMatrix[i][j] == 1 ? 0 : 255);
				fwrite(&c, 1, sizeof(unsigned char), fp);
			}
		}
		// 关闭文件
		fclose(fp);
	}

	/**
	 * @brief Get the Random Number object
	 * 
	 * @param l 随机数的左区间
	 * @param r 随机数的左区间
	 * @return int 随机数
	 */
	int getRandomNumber(int l, int r) {
		if (l > r) swap(l, r);
		return l + rand() % (r - l + 1);
	}

	/**
	 * @brief Get the Distance object
	 * 
	 * @param t1 点 1
	 * @param t2 点 2
	 * @return double 返回距离
	 */
	double getDistance(const Point &t1, const Point &t2) {
		return sqrt((t1.x - t2.x) * (t1.x - t2.x) + (t1.y - t2.y) * (t1.y - t2.y));
	}

	/**
	 * @brief Get the Distance object
	 * 
	 * @param x1 点 1 的 x 坐标
	 * @param y1 点 1 的 y 坐标
	 * @param x2 点 2 的 x 坐标
	 * @param y2 点 2 的 y 坐标
	 * @return double 返回距离
	 */
	double getDistance(double x1, double y1, double x2, double y2) {
		return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
	}

	/**
	 * @brief 在所有候选点集选择两个点
	 * 
	 * @param P 存放两个候选点
	 * @return true 在规定次数内选择到了两个候选点
	 * @return false 选择失败
	 */
	bool findTwoPoints(vector<Point> &P) {
		// 寻找次数
		int round = 0;
		while (P.size() < 2 && round++ < MAX_FIND_TWO_POINTS_ROUNDS) {
			// 产生一个在图像宽高范围内的随机点
			Point tmp(getRandomNumber(0, height), getRandomNumber(0, width));
			// 二分查找到比随机点大的候选点
			set<Point>::iterator it = candidateSet.lower_bound(tmp);
			// 查找到候选点
			if (it != candidateSet.end()) {
				P.push_back(*it);
			}
			// 已经找到两个候选点
			if (P.size() == 2) {
				// 两个候选点距离小于设置的最小两点距离
				if (getDistance(P[0], P[1]) < MIN_TWO_POINTS_DISTANCE){
					// 清空这两个候选点并重新选择
					P.clear();
				} 
				else return true;
			}
		}
		return false;
	}

	/**
	 * @brief 计算两点斜率
	 * 
	 * @param k 计算保存的斜率
	 * @param t1 点 1
	 * @param t2 点 2
	 * @return true 斜率存在
	 * @return false 斜率不存在
	 */
	bool calSlope(double &k, const Point& t1, const Point& t2) {
		if (t2.x == t1.x) return false;
		k = 1.0 * (t2.y - t1.y) / (t2.x - t1.x);
		return true;
	}

	/**
	 * @brief 计算两点中点
	 * 
	 * @param t1 点 1 
	 * @param t2 点 2
	 * @return Point 中点
	 */
	Point calMiddlePoint(const Point& t1, const Point& t2) {
		return Point((t1.x + t2.x) >> 1, (t1.y + t2.y) >> 1);
	}

	/**
	 * @brief 判断是否是孤立点（ 3 x 3 的邻域判断 ）
	 * 
	 * @param t 待判断点
	 * @return true 是孤立点
	 * @return false 不是孤立点
	 */
	bool ioslatedPoint(const Point& t) {
		for (int i = t.x - 1; i <= t.x + 1; i++) {
			for (int j = t.y - 1; j <= t.y + 1; j++) {
				if (i == t.x && j == t.y) continue;
				if (inMatrix(i, j) && originalMatrix[i][j]) return false;
			}
		}
		return true;
	}

	/**
	 * @brief 判断是否在图像矩阵内
	 * 
	 * @param x 
	 * @param y 
	 * @return true 在图像矩阵内
	 * @return false 不在图像矩阵内
	 */
	bool inMatrix(int x, int y) {
		return 0 <= x && x < height && 0 <= y && y < width;
	}

	/**
	 * @brief 判断 4 x 4 邻域是否是 半连续结构
	 * 半连续结构: 邻域边缘无像素点，中心有 2-4 个点
	 * @param lx 
	 * @param rx 
	 * @param ly 
	 * @param ry 
	 * @return true 
	 * @return false 
	 */
	bool semiContinuousStructure(int lx, int rx, int ly, int ry) {
		// 矩阵 (i, j) 有像素值为 1 的点的总数
		int sum = 0;
		// 中心 2 x 2 邻域有像素值为 1 的点的数量
		int middleSum = 0;
		for (int i = lx; i <= rx; i++) {
			for (int j = ly; j <= ry; j++) {
				if (inMatrix(i, j)) {
					sum += originalMatrix[i][j];
				}
			}
		}
		for (int i = lx + 1; i <= rx - 1; i++) {
			for (int j = ly + 1; j <= ry - 1; j++) {
				if (inMatrix(i, j)) {
					middleSum += originalMatrix[i][j];
				}
			}
		}
		return (sum - middleSum == 0) && 2 <= middleSum && middleSum <= 4;
	}

	/**
	 * @brief 判断点是否为半连续点
	 * 半连续点: 4 x 4 的邻域边缘无像素点，中心有 2-4 个点
	 * 由于是 4 x 4 的邻域, 所以四个方向的邻域只要满足一个即可
	 * @param t 
	 * @return true 
	 * @return false 
	 */
	bool semiContinuousPoint(const Point& t) {
		return semiContinuousStructure(t.x - 2, t.x + 1, t.y - 2, t.y + 1) 
			|| semiContinuousStructure(t.x - 2, t.x + 1, t.y - 1, t.y + 2) 
			|| semiContinuousStructure(t.x - 1, t.x + 2, t.y - 2, t.y + 1) 
			|| semiContinuousStructure(t.x - 1, t.x + 2, t.y - 1, t.y + 2);
	}

	/**
	 * @brief 线性拟合（ 利用 3 x 3 邻域拟合 ）
	 * 得到直线 y = ax + b
	 * @param a 
	 * @param b 
	 * @param point 中心点
	 * @return true 拟合成功
	 * @return false 拟合失败
	 */
	bool lineFitting(double& a, double& b, const Point& point) {
		double v[2][2], u[2];
		memset(v, 0, sizeof(v));
		memset(u, 0, sizeof(u));
		// 计算二元一次方程系数
		v[0][0] = 9;
		for (int i = point.x - 1; i <= point.x + 1; i++) {
			for (int j = point.y - 1; j <= point.y + 1; j++) {
				if (inMatrix(i, j)) {
					v[0][1] += j * originalMatrix[i][j];
					v[1][0] += j * originalMatrix[i][j];
					v[1][1] += j * j * originalMatrix[i][j];
				}
			}
		}
		for (int i = point.x - 1; i <= point.x + 1; i++) {
			for (int j = point.y - 1; j <= point.y + 1; j++) {
				if (inMatrix(i, j)) {
					u[0] += i * originalMatrix[i][j];
					u[1] += i * j * originalMatrix[i][j];
				}
			}
		}
		// 二元一次方程无解
		if (v[0][1] * v[1][0] - v[0][0] * v[1][1] == 0) return false;
		// 二元一次方程有解
		a = (u[0] * v[1][0] - u[1] * v[0][0]) / (v[0][1] * v[1][0] - v[0][0] * v[1][1]);
		b = (u[1] * v[0][1] - u[0] * v[1][1]) / (v[0][1] * v[1][0] - v[0][0] * v[1][1]);
		return true;
	}

	/**
	 * @brief 计算两斜率直线的角度
	 * 
	 * @param k1 
	 * @param k2 
	 * @return double 角度（ 0 - 360 ）
	 */
	double calAngle(double k1, double k2) {
		// 当斜率不存在时, 返回 double 的最大值
		if (1.0 + k1 * k2 == 0) return DBL_MAX;
		return atan(abs(k1 - k2) / (1.0 + k1 * k2)) * 180.0 / PI;
	}

	/**
	 * @brief 判断候选三点是否共圆
	 * 
	 * @param k 拟合直线的斜率
	 * @param point 第三个候选点位置
	 * @param M 第一第二个点的中点
	 * @return true 
	 * @return false 
	 */
	bool judgeDirection(int k, const Point& point, const Point& M) {
		// count = 直线上的点数 + 同侧点数 - 异侧点数
		int count = 0;
		for (int i = point.x - 1; i <= point.x + 1; i++) {
			for (int j = point.y - 1; j <= point.y + 1; j++) {
				if (inMatrix(i, j) && originalMatrix[i][j] == 1) {
					double status = (j - point.y - k * (i - point.x)) * (M.y - point.y - k * (M.x - point.x));
					count += (status >= 0 ? 1 : -1);
				}
			}
		}
		return count > JUDGE_DIRECTION_THRESHOLD;
	}

	/**
	 * @brief 寻找第三个候选点位置
	 * 
	 * @param P 存放三个候选点
	 * @param M 第一第二候选点中点位置
	 * @param k 第一二候选点的斜率
	 * @return true 寻找到第三点
	 * @return false 
	 */
	bool findThirdPoint(vector<Point>& P, const Point& M, double& k) {
		// Y - M.y = -1.0 / k * (X - M.x)
		// x 轴正方向寻找
		//printf("  {%d, %d}, {%d, %d}\n", P[0].x, P[0].y, P[1].x, P[1].y);
		for (int x = M.x + 1; x < height; x++) {
			int y = -1.0 / k * (x - M.x) + M.y;
			Point point(x, y);
			// 在图像矩阵范围内 && 该位置存在有效点 && 不是孤立点 && 不是半连续点
			if (0 <= y && y < width && originalMatrix[x][y] && !ioslatedPoint(point) && !semiContinuousPoint(point)) {
				double a, b;
				// 直线拟合
				if (!lineFitting(a, b, point)) continue;
				// 计算角度误差是否满足条件
				if (calAngle(k, a) > ANGLE_ERROR) continue;
				// 判断三点共圆
				if (!judgeDirection(a, point, M)) continue;
				// 将第三点加入 P
				P.push_back(point);
				return true;
			}
		}
		// x 轴负方向寻找
		for (int x = M.x - 1; x >= 0; x--) {
			int y = -1.0 / k * (x - M.x) + M.y;
			Point point(x, y);
			// 在图像矩阵范围内 && 该位置存在有效点 && 不是孤立点 && 不是半连续点
			if (0 <= y && y < width && originalMatrix[x][y] && !ioslatedPoint(point) && !semiContinuousPoint(point)) {
				double a, b;
				// 直线拟合
				if (!lineFitting(a, b, point)) continue;
				// 计算角度误差是否满足条件
				if (calAngle(k, a) > ANGLE_ERROR) continue;
				// 判断三点共圆
				if (!judgeDirection(a, point, M)) continue;
				// 将第三点加入 P
				P.push_back(point);
				return true;
			}
		}
		return false;
	}

	/**
	 * @brief 通过三点计算圆的参数
	 * 
	 * @param a 圆心坐标 x
	 * @param b 圆心坐标 y
	 * @param r 圆半径
	 * @param P 三个候选点
	 * @return true 三个候选点不共线可以计算出圆参数
	 * @return false 三个候选点共线
	 */
	bool calCandidateCycleParameter(double& a, double& b, double& r, vector<Point>& P) {
		double x1 = P[0].x, y1 = P[0].y;
		double x2 = P[1].x, y2 = P[1].y;
		double x3 = P[2].x, y3 = P[2].y;
		double A = x1 * (y2 - y3) - y1 * (x2 - x3) + x2 * y3 - x3 * y2;
		double B = (x1 * x1 + y1 * y1) * (y3 - y2) + (x2 * x2 + y2 * y2) * (y1 - y3) + (x3 * x3 + y3 * y3) * (y2 - y1);
		double C = (x1 * x1 + y1 * y1) * (x2 - x3) + (x2 * x2 + y2 * y2) * (x3 - x1) + (x3 * x3 + y3 * y3) * (x1 - x2);
		double D = (x1 * x1 + y1 * y1) * (x3 * y2 - x2 * y3) + (x2 * x2 + y2 * y2) * (x1 * y3 - x3 * y1) + (x3 * x3 + y3 * y3) * (x2 * y1 - x1 * y2);
		if (A == 0) return false;
		a = -B / (2.0 * A);
		b = -C / (2.0 * A);
		r = sqrt((B * B + C * C - 4.0 * A * D) / (4.0 * A * A));
		return true;
	}

	/**
	 * @brief 得到大正方形和小正方形之间的有效点集
	 * 
	 * @param candidatePoints 记录有效点集
	 * @param a 圆心 x 坐标
	 * @param b 圆心 y 坐标
	 * @param r1 大正方形边长
	 * @param r2 小正方形边长
	 */
	void getCandidatePoints(vector<Point>& candidatePoints, double a, double b, double r1, double r2) {
		int outlx = a - r1 / 2.0;
		int outrx = a + r1 / 2.0;
		int outly = b - r1 / 2.0;
		int outry = b + r1 / 2.0;
		int inlx = a - r2 / 2.0;
		int inrx = a + r2 / 2.0;
		int inly = b - r2 / 2.0;
		int inry = b + r2 / 2.0;
		//printf("	out(lx:%d, rx:%d, ly:%d, ry:%d)\n", outlx, outrx, outly, outry);
		//printf("	in (lx:%d, rx:%d, ly:%d, ry:%d)\n", inlx, inrx, inly, inry);
		for (int i = outlx; i <= outrx; i++) {
			for (int j = outly; j <= outry; j++) {
				if (inlx <= i && i <= inrx && inly <= j && j <= inry) {
					j = inry;
					continue;
				}
				if (inMatrix(i, j) && originalMatrix[i][j]) {
					candidatePoints.push_back(Point(i, j));
				}
			}
		}
	}

	/**
	 * @brief 判断是否是真圆
	 * 判断依据: 图像矩阵中在圆上的点的数量超过一定比例
	 * @param a 圆心 x 坐标
	 * @param b 圆心 y 坐标
	 * @param r 圆半径
	 * @param realCyclePoints 图像矩阵中真圆点集
	 * @return true 
	 * @return false 
	 */
	bool judgeRealCycle(double a, double b, double r, vector<Point>& realCyclePoints) {
		// 获取大小正方形之间的有效点集
		vector<Point> candidatePoints;
		getCandidatePoints(candidatePoints, a, b, 2.0 * r + 2, 1.1 * r);
		// 有效点集小于圆上所需点数
		if (candidatePoints.size() < PROPORTIONALITY_FACTOR * 2.0 * PI * r) return false;
		// 遍历有效点集是否在圆上并保存至真圆点集
		for (auto& point : candidatePoints) {
			// 计算有效点和圆心距离
			double dis = getDistance(a, b, point.x, point.y);
			// 距离小于一定阈值
			if (abs(getDistance(a, b, point.x, point.y) - r) < DISTANCE_ERROR) {
				realCyclePoints.push_back(point);
			}
		}
		// 真圆点数大于一定比例
		return realCyclePoints.size() >= PROPORTIONALITY_FACTOR * 2.0 * PI * r;
	}

	/**
	 * @brief 记录真圆参数
	 * 
	 * @param a 圆心 x 坐标
	 * @param b 圆心 y 坐标
	 * @param r 圆半径
	 */
    void recordRealCycleParameter(double a, double b, double r){
        resultCycleSet.push_back(Cycle(a, b, r));
    }

	/**
	 * @brief 清除图像矩阵真圆点集
	 * 
	 * @param realCyclePoints 
	 */
	void clearOriginalMatrixRealCyclePoints(const vector<Point>& realCyclePoints) {
		for (auto& point : realCyclePoints) {
			originalMatrix[point.x][point.y] = 0;
			candidateSet.erase(point);
		}
	}

	/**
	 * @brief 将真圆集合映射到最终矩阵
	 * 
	 */
    void resultCycleSetToFinalMatrix(){
        for (auto& cycle : resultCycleSet) {
            for (double angle = 0; angle < 360; angle += 0.1) {
                int x = cycle.x + cycle.r * cos(angle * PI / 180.0);
                int y = cycle.y + cycle.r * sin(angle * PI / 180.0);
                finalMatrix[x][y] = 1;
            }
        }
    }

	/**
	 * @brief 运行我的多圆检测算法
	 * 
	 */
	void run() {
		srand(time(0));
		int round = 0;
		while (round++ < MAX_ROUNDS) {
			//printf("round = %d ========================= \n", round);
			vector<Point> P;
			// 寻找两点
			if (!findTwoPoints(P)) continue;
			// 计算斜率
			double k;
			if (!calSlope(k, P[0], P[1])) continue;
			// 计算中点
			Point M = calMiddlePoint(P[0], P[1]);
			// 寻找第三点
			if (!findThirdPoint(P, M, k)) continue;
			// 计算圆参数
			double a, b, r;
			if (!calCandidateCycleParameter(a, b, r, P)) continue;
			// 判断是否为真圆
			vector<Point> realCyclePoints;
			if (!judgeRealCycle(a, b, r, realCyclePoints)) continue;
            // 记录真圆圆参数
            recordRealCycleParameter(a, b, r);
			// 清除原图真圆点集
			clearOriginalMatrixRealCyclePoints(realCyclePoints);
			round = 0;
		}
		// 真圆集映射到结果矩阵
        resultCycleSetToFinalMatrix();
	}
};