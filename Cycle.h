#pragma once

/**
 * @brief 圆类
 * 
 */
class Cycle {
public:
	// 圆心 x, y
	double x, y;
	// 圆半径
	double r;

	/**
	 * @brief Construct a new Cycle object
	 * 
	 * @param x 
	 * @param y 
	 * @param r 
	 */
	Cycle(double x, double y, double r) : x(x), y(y), r(r) {}

	/**
	 * @brief 打印圆参数
	 * 
	 */
	void show() {
		printf("{%f, %f, %f}\n", x, y, r);
	}
};