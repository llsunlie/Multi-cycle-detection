#pragma once

/**
 * @brief 二维平面点类
 * 
 */
class Point {
public:
	int x;
	int y;
	
	/**
	 * @brief Construct a new Point object
	 * 
	 * @param x 
	 * @param y 
	 */
	Point(int x, int y) : x(x), y(y) {}
	
	/**
	 * @brief 重载小于运算符 x 越小优先级越高, y 越小优先级越高
	 * 
	 * @param another 
	 * @return true 
	 * @return false 
	 */
	bool operator < (const Point& another) const {
		if (this->x == another.x) return this->y < another.y;
		return this->x < another.x;
	}
};