/*
 * CalcDistanceOfFoldLine.cpp
 *
 *  Created on: 2015年9月1日
 *      Author: netbeen
 */


/*
 * ref1: http://blog.csdn.net/angelazy/article/details/38489293	点到线段的最短距离
 * ref2: http://www.docin.com/p-476754366.html	空间中两条线段之间的最短距离.doc
 */

#include <vector>
#include <assert.h>
#include <cmath>
#include <iostream>
#include <climits>

////3维点
class Point3f {
public:
	const double x;
	const double y;
	const double z;
	Point3f(double x, double y, double z) :
			x(x), y(y), z(z) {
	}
};

//3维向量
class Vector3f {
public:
	const double x;
	const double y;
	const double z;
	Vector3f(double x, double y, double z) :
			x(x), y(y), z(z) {
	}
};

class CalcDistanceOfFoldLine {
private:
	std::vector<Point3f> foldLine1;
	std::vector<Point3f> foldLine2;

	//采样函数
	std::vector<Point3f> inline sampling(Point3f point1, Point3f point2, int num = 100) {
		std::vector<Point3f> result;

		float stepX = (point2.x - point1.x) / num;
		float stepY = (point2.y - point1.y) / num;
		float stepZ = (point2.z - point1.z) / num;

		for (int i = 0; i < num; i++) {
			Point3f point = Point3f(point1.x + stepX * i, point1.y + stepY * i, point1.z + stepZ * i);
			result.push_back(point);
		}

		return result;
	}

	//求两个点的距离
	float inline getDistance(Point3f point1, Point3f point2) {
		return sqrt(pow(point1.x - point2.x, 2) + pow(point1.y - point2.y, 2) + pow(point1.z - point2.z, 2));
	}

	//求两个向量的点积
	double  scalarProduct3f(double a1, double a2, double a3, double b1, double b2, double b3){
		return a1*b1+a2*b2+a3*b3;
	}

	//求向量的模
	double magnitude(double a1, double a2, double a3){
		return sqrt(pow(a1,2)+pow(a2,2)+pow(a3,2));
	}

	//求两个点的距离
	double distanceOfPointToPoint(Point3f point1, Point3f point2){
		return sqrt(pow(point1.x-point2.x,2)+pow(point1.y-point2.y,2)+pow(point1.z-point2.z,2));
	}

	//求点到线段的距离
	double distanceOfPointToLine(Point3f point, Point3f lineStart, Point3f lineEnd){
		Vector3f ab = Vector3f(lineEnd.x-lineStart.x,lineEnd.y-lineStart.y,lineEnd.z-lineStart.z);
		Vector3f ap = Vector3f(point.x-lineStart.x,point.y-lineStart.y,point.z-lineStart.z);
		double r =this->scalarProduct3f(ab.x,ab.y,ab.z,ap.x,ap.y,ap.z) / this->magnitude(ab.x,ab.y,ab.z);

		if(r>1){
			return this->distanceOfPointToPoint(point, lineEnd);
		}else if(r < 0){
			return this->distanceOfPointToPoint(point, lineStart);
		}else{
			double Mab = this->magnitude(ab.x,ab.y,ab.z);	//ab的模
			double Map = this->magnitude(ap.x,ap.y,ap.z);	//ap的模
			return sqrt(pow(Map,2)-pow(Mab*r,2));
		}
	}

public:
	//输入两个折线段集，返回一个五元素vector（距离、第一个折线段首部index，第一个折线段尾部部index，第二个折线段首部index，第二个折线段尾部index）
	std::vector<double> calc(std::vector<Point3f> foldLine1, std::vector<Point3f> foldLine2) {
		assert(foldLine1.size() >= 2);
		assert(foldLine2.size() >= 2);

		std::vector<double> result;
		float minDistance = INT_MAX;
		int startIndexI = -1;
		int endIndexI = -1;
		int startIndexJ = -1;
		int endIndexJ = -1;

		for (int i = 0; i < foldLine1.size() - 1; i++) {
			//std::vector<Point3f> currentSampleI = this->sampling(foldLine1.at(i), foldLine1.at(i + 1));
			for (int j = 0; j < foldLine2.size() - 1; j++) {
				double x1 = foldLine1.at(i).x;
				double y1 = foldLine1.at(i).y;
				double z1 = foldLine1.at(i).z;
				double x2 = foldLine1.at(i + 1).x;
				double y2 = foldLine1.at(i + 1).y;
				double z2 = foldLine1.at(i + 1).z;
				double x3 = foldLine2.at(i).x;
				double y3 = foldLine2.at(i).y;
				double z3 = foldLine2.at(i).z;
				double x4 = foldLine2.at(i + 1).x;
				double y4 = foldLine2.at(i + 1).y;
				double z4 = foldLine2.at(i + 1).z;

				//std::vector<Point3f> currentSampleJ = this->sampling(foldLine2.at(j), foldLine2.at(j + 1));
				/*for (Point3f elemI : currentSampleI) {
				 for (Point3f elemJ : currentSampleJ) {
				 float distance = this->getDistance(elemI,elemJ);
				 if(distance < minDistance){
				 minDistance = distance;
				 startIndexI = i;
				 endIndexI = i+1;
				 startIndexJ = j;
				 endIndexJ = j+1;
				 }
				 }
				 }*/

				double a1 = pow(x2 - x1, 2) + pow(y2 - y1, 2) + pow(z2 - z1, 2);
				double b1 = -((x2 - x1) * (x4 - x3) + (y2 - y1) * (y4 - y3) + (z2 - z1) * (z4 - z3));
				double c1 = (x1 - x2) * (x1 - x3) + (y1 - y2) * (y1 - y3) + (z1 - z2) * (z1 - z3);
				double a2 = -((x2 - x1) * (x4 - x3) + (y2 - y1) * (y4 - y3) + (z2 - z1) * (z4 - z3));
				double b2 = pow(x4 - x3, 2) + pow(y4 - y3, 2) + pow(z4 - z3, 2);
				double c2 = (x1 - x3) * (x4 - x3) + (y1 - y3) * (y4 - y3) + (z1 - z3) * (z4 - z3);
				double t = (a1 * c2 - a2 * c1) / (a1 * b2 - a2 * b1);
				double s = (c1 - b1 * t) / a1;

				//std::cout << "s=" << s << " t=" << t << std::endl;
				double localDistance;
				if (0 <= s && s <= 1 && 0 <= t && t <= 1) {
					double X = x1 + s * (x2 - x1);
					double Y = y1 + s * (y2 - y1);
					double Z = z1 + s * (z2 - z1);

					double U = x3 + s * (x4 - x3);
					double V = y3 + s * (y4 - y3);
					double W = z3 + s * (z4 - z3);
					localDistance = sqrt(pow(X - U, 2) + pow(Y - V, 2) + pow(Z - W, 2));
					//std::cout << "localDistance: " << localDistance << std::endl;
				}else{
					double distanceAToCD = this->distanceOfPointToLine(foldLine1.at(i),foldLine2.at(j),foldLine2.at(j+1));
					//std::cout << "distanceAToCD: " << distanceAToCD << std::endl;
					double distanceBToCD = this->distanceOfPointToLine(foldLine1.at(i+1),foldLine2.at(j),foldLine2.at(j+1));
					//std::cout << "distanceBToCD: " << distanceBToCD << std::endl;
					double distanceCToAB = this->distanceOfPointToLine(foldLine2.at(j),foldLine1.at(i),foldLine1.at(i+1));
					//std::cout << "distanceCToAB: " << distanceCToAB << std::endl;
					double distanceDToAB = this->distanceOfPointToLine(foldLine2.at(j+1),foldLine1.at(i),foldLine1.at(i+1));
					//std::cout << "distanceDToAB: " << distanceDToAB << std::endl;

					localDistance = std::min(std::min(distanceAToCD,distanceBToCD),std::min(distanceCToAB,distanceDToAB));
					//std::cout << "localDistance: " << localDistance << std::endl;
				}

				if(localDistance < minDistance){
					minDistance = localDistance;
					 startIndexI = i;
					 endIndexI = i+1;
					 startIndexJ = j;
					 endIndexJ = j+1;
				}

			}
		}

		result.push_back(minDistance);
		result.push_back(startIndexI);
		result.push_back(endIndexI);
		result.push_back(startIndexJ);
		result.push_back(endIndexJ);
		return result;
	}

};

int main() {

	std::vector<Point3f> foldLine1 { Point3f(100, 100, 100),Point3f(1, 1, 1), Point3f(0, 0, 1) };
	std::vector<Point3f> foldLine2 { Point3f(-100, -100, -100),Point3f(0, 1, 0), Point3f(1, 0, 0) };

	CalcDistanceOfFoldLine s;
	std::vector<double> result = s.calc(foldLine1, foldLine2);

	std::cout << "Distance:" << std::endl;
	std::cout << result.front() << std::endl;
	std::cout << "Index:" << std::endl;
	std::cout << result.at(1) << std::endl;
	std::cout << result.at(2) << std::endl;
	std::cout << result.at(3) << std::endl;
	std::cout << result.at(4) << std::endl;
	return 0;
}
