//----------------------------------------------------------
/// @file
/// @brief    HPCAnswer.hpp の実装 (解答記述用ファイル)
/// @author   ハル研究所プログラミングコンテスト実行委員会
///
/// @copyright  Copyright (c) 2013 HAL Laboratory, Inc.
/// @attention  このファイルの利用は、同梱のREADMEにある
///             利用条件に従ってください

//----------------------------------------------------------

// Answer.cpp 専用のインクルードファイルです。
// 別のファイルをインクルードした場合、評価時には削除されます。
#include "HPCAnswerInclude.hpp"
#include <iostream>
float ADD(float a, float b){
	const double EPS = 1e-3;
	if (hpc::Math::Abs(a + b) <= EPS*(hpc::Math::Abs(a) + hpc::Math::Abs(b))) return 0;
	return a + b;
}
double ABS(double a){
	if (a >= 0) return a; return -a;
}
double ADD(double a, double b){
	const double EPS = 1e-10;
	if (ABS(a + b) <= EPS*(ABS(a) + ABS(b))) return 0;
	return a + b;
}
struct pair{
	float x; int y;
	pair() {}
	pair(float x, int y) : x(x), y(y){}
};
struct pair_ii{
	int x, y;
	pair_ii() {}
	pair_ii(int x, int y) : x(x), y(y){}
};
struct P{
	float x,y;
	P() {}
	P(float x, float y) : x(x), y(y) {}
	P operator+(P p) { return P(ADD(x, p.x), ADD(y, p.y)); }
	P operator-(P p) { return P(ADD(x, -p.x), ADD(y, -p.y)); }
	P operator*(float d){ return P(x*d, y*d); }
	float dot(P p){ return ADD(x*p.x, y*p.y); }
	float det(P p){ return ADD(x*p.y, -y*p.x); }
};
class VEC{
		/*
	double add(double a, double b){
		const double EPS = 1e-10;
		if (abs(a + b) <= EPS*(abs(a) + abs(b))) return 0;
		return a + b;
	}
	*/
public:
	double x, y;
	VEC() { x = 0, y = 0; };
	VEC(double a, double b){ x = a, y = b; };
	VEC operator+(VEC p){ return VEC(ADD(x, p.x), ADD(y, p.y)); }
	VEC operator-(VEC p){ return VEC(ADD(x, -p.x), ADD(y, -p.y)); }
	VEC operator*(double d){ return VEC(x*d, y*d); }
	VEC operator/(double d){ return VEC(x / d, y / d); }
	double dot(VEC p){ return ADD(x*p.x, y*p.y); }
	double det(VEC p){ return ADD(x*p.y, -y*p.x); }
	double squareLength(){ return ADD(x*x, y*y); }
	/*
	double Length(){ return sqrt(ADD(x*x, y*y)); }
	double rotSign(VEC aTarget){
		return atan2(x * aTarget.y - y * aTarget.x, dot(aTarget));
	}
	double cos(VEC v){
		return dot(v) / Length() * v.Length();
	}
	*/
};
/// プロコン問題環境を表します。
namespace hpc {

	//----------------------------------------------------------
	/// 各ステージ開始時に呼び出されます。
	///
	/// この関数を実装することで、各ステージに対して初期設定を行うことができます。
	///
	/// @param[in] aStage 現在のステージ。
	void Answer::Init(const Stage& aStage)
	{
	}
	const int HoleNum(const Vec2& v, const HoleCollection& aHoles)
	{
		int index = 0;
		for (; index < aHoles.count(); index++){
			if (aHoles[index].isIn(v)) break;
		}
		if (index == aHoles.count()) return -1;
		else return index;
	}
	void qsort_P(pair* a, int left, int right){
		int pl = left;
		int pr = right;
		float x = a[(pl + pr) / 2].x;
		while (pl <= pr){
			while (a[pl].x < x) pl++;
			while (a[pr].x > x) pr--;
			if (pl <= pr){
				pair tmp = a[pl]; a[pl] = a[pr], a[pr] = tmp;
				pl++, pr--;
			}
		}
		if (left < pr) qsort_P(a, left, pr);
		if (pl < right) qsort_P(a, pl, right);
	}
	void qsort_P(pair_ii* a, int left, int right){
		int pl = left;
		int pr = right;
		int x = a[(pl + pr) / 2].x;
		while (pl <= pr){
			while (a[pl].x < x) pl++;
			while (a[pr].x > x) pr--;
			if (pl <= pr){
				pair_ii tmp = a[pl]; a[pl] = a[pr], a[pr] = tmp;
				pl++, pr--;
			}
		}
		if (left < pr) qsort_P(a, left, pr);
		if (pl < right) qsort_P(a, pl, right);
	}
	void qsort_sp(int* a, int left, int right, const Stage& aStage){
		int pl = left;
		int pr = right;
		float x = aStage.items()[a[(pl + pr) / 2]].pos().x;
		float y = aStage.items()[a[(pl + pr) / 2]].pos().y;
		while (pl <= pr){
			while (ADD(aStage.items()[a[pl]].pos().x, -x) < 0) pl++;
			while (ADD(aStage.items()[a[pl]].pos().x, -x) == 0 && ADD(aStage.items()[a[pl]].pos().y, -y)<0) pl++;
			while (ADD(aStage.items()[a[pr]].pos().x, -x) > 0) pr--;
			while (ADD(aStage.items()[a[pr]].pos().x, -x) == 0 && ADD(aStage.items()[a[pr]].pos().y, -y)>0) pr--;
			if (pl <= pr){
				int tmp = a[pl]; a[pl] = a[pr], a[pr] = tmp;
				pl++, pr--;
			}
		}
		if (left < pr) qsort_sp(a, left, pr, aStage);
		if (pl < right) qsort_sp(a, pl, right, aStage);
	}
	//xorShift
	unsigned long xor128(){
		static unsigned long x = 123456789, y = 362436069;
		static unsigned long z = 521288629, w = 88675123;
		unsigned long t;
		t = (x ^ (x << 11)); x = y, y = z, z = w;
		return (w = (w ^ (w >> 19)) ^ (t ^ (t >> 8)));
	}
	bool IsHit(const Circle& aC0, const Circle& aC1)
	{
		const float squareDist = aC0.pos().squareDist(aC1.pos());
		return squareDist <= (aC0.radius() + aC1.radius()) * (aC0.radius() + aC1.radius());
	}
	bool WillHit(const Circle& aC0, const Circle& aC1, const Vec2& aP1)
	{
		// 動いていなければ静止している円の判定
		if (aC1.pos() == aP1) {
			return IsHit(aC0, aC1);
		}

		// 実質的に半径 aC0.radius + aC1.radius の円と線分 [aC1.pos, aP1.pos] の衝突判定と
		// 同等になる。
		const Circle c(aC0.pos(), aC0.radius() + aC1.radius());
		const Vec2 seg = aP1 - aC1.pos();
		const Vec2 c1ToC0 = c.pos() - aC1.pos();
		const float dist = Math::Abs(seg.cross(c1ToC0)) / seg.length();
		// 距離が c.radius より遠ければ衝突しない
		if (dist > c.radius()) {
			return false;
		}
		// 線分の延長線上で交差していないか調べる。
		const Vec2 p1ToC0 = c.pos() - aP1;
		// それぞれの点が円の反対方向にあれば衝突
		if (c1ToC0.dot(seg) * p1ToC0.dot(seg) <= 0.0f) {
			return true;
		}
		// 半径が大きければ衝突
		if (c.radius() >= c1ToC0.length() || c.radius() >= p1ToC0.length()) {
			return true;
		}
		return false;
	}
	bool on_seg(Vec2 s1, Vec2 f1, Vec2 q){
		return ADD((s1 - q).x*(f1 - q).y, -(s1 - q).y*(f1 - q).x) == 0 && ADD((s1 - q).x*(f1 - q).x, (s1 - q).y*(f1 - q).y) <= 0;
	}
	bool on_seg(VEC s1, VEC f1, VEC q){
		return ADD((s1 - q).x*(f1 - q).y, -(s1 - q).y*(f1 - q).x) == 0 && ADD((s1 - q).x*(f1 - q).x, (s1 - q).y*(f1 - q).y) <= 0;
	}
	//線分(s1-f1)と線分(s2-f2)の交点
	Vec2 intersection(Vec2  s1, Vec2 f1, Vec2 s2, Vec2 f2){
		return s1 + (f1 - s1)*((f2 - s2).cross(s2 - s1) / (f2 - s2).cross(f1 - s1));
	}
	VEC intersection(VEC  s1, VEC f1, VEC s2, VEC f2){
		return s1 + (f1 - s1)*((f2 - s2).det(s2 - s1) / (f2 - s2).det(f1 - s1));
	}
	bool isCross(Vec2 s1, Vec2 f1, Vec2 s2, Vec2 f2){
		if (ADD((f1 - s1).x*(f2 - s2).y, -(f1 - s1).y*(f2 - s2).x) == 0){
			return false;
			return on_seg(s1, f1, s2) || on_seg(s1, f1, f2) || on_seg(s2, f2, s1) || on_seg(s2, f2, f1);
		}
		else{
			Vec2 r = intersection(s1, f1, s2, f2);
			return on_seg(s1, f1, r) && on_seg(s2, f2, r);
		}
	}
	bool isCross(VEC s1, VEC f1, VEC s2, VEC f2){
		if (ADD((f1 - s1).x*(f2 - s2).y, -(f1 - s1).y*(f2 - s2).x) == 0){
			return false;
			return on_seg(s1, f1, s2) || on_seg(s1, f1, f2) || on_seg(s2, f2, s1) || on_seg(s2, f2, f1);
		}
		else{
			VEC r = intersection(s1, f1, s2, f2);
			return on_seg(s1, f1, r) && on_seg(s2, f2, r);
		}
	}
	//O(n)
	const int FindNearestItemNumGreedy(const ItemCollection& items, const Stage& aStage)
	{
		int MIN = 100000; int tmp = 0;
		//未取得 最も近い
		for (int index = 0; index < items.count(); ++index) {
			if (!items[index].isRemoved()){
				float DIST = (items[index].pos() - aStage.player().pos()).length();
				Vec2 playerToItem = aStage.items()[index].pos() - aStage.player().pos();
				//正式にはMath::Ceilをつけるべき、ただしスコアは変わらない
				int turn = (DIST / Parameter::PLAYER_SPEED);
				int rotTurn = Math::Abs(aStage.player().vec().rotSign(playerToItem)) / Parameter::PLAYER_ANGLE_RATE;
				for (int j = rotTurn; j >= 0; j--){
					const Circle Item(aStage.items()[index].pos(), aStage.items()[index].radius());
					const Circle PLAYER(aStage.player().pos(), aStage.player().region().radius());
					//ActiveCollisionはそのまま進んだときにItemを取れるかの判定
					float arg = aStage.player().arg();
					if (aStage.player().vec().rotSign(playerToItem) >= 0)arg += rotTurn * Parameter::PLAYER_ANGLE_RATE;
					else arg -= rotTurn * Parameter::PLAYER_ANGLE_RATE;
					Vec2 playerVec = Vec2(Math::Cos(arg), Math::Sin(arg));
					const Vec2 Velocity = Parameter::PLAYER_SPEED * playerVec;
					const bool ActiveCollision = WillHit(Item, PLAYER, aStage.player().pos() + turn * Velocity);
					if (!ActiveCollision) break;
					rotTurn--;
				}
				turn += rotTurn;
				if (MIN > turn) MIN = turn, tmp = index;
			}
		}
		if (MIN != 100000) return tmp; else return -1;
	}

	const int FindNearestItemNumGreedyDist(const ItemCollection& items, const Stage& aStage)
	{
		float MIN = 100000; int tmp = 0;
		//未取得 最も近い
		for (int index = 0; index < items.count(); ++index) {
			if (!items[index].isRemoved()){
				float DIST = (items[index].pos() - aStage.player().pos()).length();
				if (MIN > DIST) MIN = DIST, tmp = index;
			}
		}
		if (MIN != 100000) return tmp; else return -1;
	}

	const int FindeNearestItemNumWithConvexHull(const ItemCollection& items, const Stage& aStage){
		//convex hull O(lg n)
		P v[110]; int v_size = 0;
		int num[110]; int num_size = 0;
		P p[110]; int p_size = 0;
		bool used[110] = { 0 };
		int NumtoItemNum[110], tmp_num[110]; int tmp_num_size = 0;
		for (int i = 0; i < items.count(); i++)	if (!items[i].isRemoved()) tmp_num[tmp_num_size++] = i;
		qsort_sp(tmp_num, 0, tmp_num_size - 1, aStage);
		for (int i = 0; i < tmp_num_size; i++){
			p[p_size] = P(items[tmp_num[i]].pos().x, items[tmp_num[i]].pos().y);
			NumtoItemNum[p_size] = tmp_num[i];
			p_size++;
		}
		//num: v[i]のi -> p[i]のi に対応させる
		//NumtoItemNum: num[i]のi -> Items[i]のi に対応させる
		int ans[110]; int ans_size = 0;
		for (int i = 0; i < p_size; i++){
			if (used[i]) continue;
			for (bool loop = true; loop;){
				loop = false;
				if (v_size < 2) { v[v_size++] = p[i], num[num_size++] = i, used[i] = true; continue; }
				if (v_size == 2){
					if ((v[1] - v[0]).det(p[i] - v[1]) <= 0) v[v_size++] = p[i], num[num_size++] = i, used[i] = true;
					else used[num[1]] = false, num_size--, v_size--, loop = true;
					continue;
				}
				int len = v_size;
				if ((v[len - 1] - v[len - 2]).det(p[i] - v[len - 1]) <= 0) v[v_size++] = p[i], num[num_size++] = i, used[i] = true;
				else used[num[len - 1]] = false, num_size--, v_size--, loop = true;
			}
		}
		for (int i = 0; i < num_size; i++) ans[ans_size++] = NumtoItemNum[num[i]];
		//
		v_size = 0, num_size = 0; for (int i = 0; i < p_size; i++) used[i] = false;
		for (int i = p_size - 1; i >= 0; i--){
			if (used[i]) continue;
			for (bool loop = true; loop;){
				loop = false;
				if (v_size < 2) { v[v_size++] = p[i], num[num_size++] = i, used[i] = true; continue; }
				if (v_size == 2){
					if ((v[1] - v[0]).det(p[i] - v[1]) <= 0) v[v_size++] = p[i], num[num_size++] = i, used[i] = true;
					else used[num[1]] = false, num_size--, v_size--, loop = true;
					continue;
				}
				int len = v_size;
				if ((v[len - 1] - v[len - 2]).det(p[i] - v[len - 1]) <= 0) v[v_size++] = p[i], num[num_size++] = i, used[i] = true;
				else used[num[len - 1]] = false, num_size--, v_size--, loop = true;
			}
		}
		for (int i = 1; i < num_size - 1; i++) ans[ans_size++] = NumtoItemNum[num[i]];
		//vに頂点が含まれている(反時計周り)
		float MIN = 100000; int MIN_NUM = -1;
		for (int i = 0; i < ans_size; i++){
			float TMP = (items[ans[i]].pos() - aStage.player().pos()).length();
			if (MIN>TMP) MIN = TMP, MIN_NUM = ans[i];
		}
		int GreedySolution = FindNearestItemNumGreedy(items, aStage);
		if (MIN_NUM == GreedySolution) return MIN_NUM;
		int cnt = 0;
		const int n = items.count();
		for (int i = 0; i < n; i++) if (aStage.items()[i].isRemoved()) used[i] = true, cnt++; else used[i] = false;
		float ConvexHullDist = MIN; int curr = MIN_NUM; used[curr] = true;
		for (int j = 0; j < n - cnt; j++){
			MIN = 100000; int NEXT = -1;
			for (int i = 0; i < n; i++){
				if (used[i]) continue;
				float DIST = (items[i].pos() - items[curr].pos()).length();
				if (MIN>DIST) MIN = DIST, NEXT = i;
			}
			if (NEXT != -1) { ConvexHullDist += MIN; curr = NEXT; used[curr] = true; }
		}

		for (int i = 0; i < n; i++) if (aStage.items()[i].isRemoved()) used[i] = true; else used[i] = false;
		curr = GreedySolution; used[curr] = true;
		float GreedyDist = (items[curr].pos() - aStage.player().pos()).length();
		for (int j = 0; j < n - cnt; j++){
			MIN = 100000; int NEXT = -1;
			for (int i = 0; i < n; i++){
				if (used[i]) continue;
				float DIST = (items[i].pos() - items[curr].pos()).length();
				if (MIN>DIST) MIN = DIST, NEXT = i;
			}
			if (NEXT != -1) { GreedyDist += MIN; curr = NEXT; used[curr] = true; }
		}
		if (GreedyDist < ConvexHullDist) return GreedySolution;
		return MIN_NUM;
	}

	Action SearchWayToEdge(const Stage& aStage, const int hole_num, const int item_num, int PredictTurn)
	{
		//目標:item_num, 穴:hole_num
		pair_ii dist_HoleEdge[8]; float MIN = 100000;
		Vec2 playerDir = aStage.player().vec(); const float CONST = 0.05;
		Vec2 d[8]; const int AppNum = 4;
		d[0] = aStage.holes()[hole_num].pointA() + CONST*Vec2(-1.0, 1.0);
		d[1] = aStage.holes()[hole_num].pointB() + CONST*Vec2(-1.0, -1.0);
		d[2] = aStage.holes()[hole_num].pointC() + CONST*Vec2(1.0, -1.0);
		d[3] = aStage.holes()[hole_num].pointD() + CONST*Vec2(1.0, 1.0);
		d[4] = (aStage.holes()[hole_num].pointA() + aStage.holes()[hole_num].pointB()) / 2 + Vec2(-1.0, 0.0) * CONST;
		d[5] = (aStage.holes()[hole_num].pointB() + aStage.holes()[hole_num].pointC()) / 2 + Vec2(0.0, -1.0) * CONST;
		d[6] = (aStage.holes()[hole_num].pointC() + aStage.holes()[hole_num].pointD()) / 2 + Vec2(1.0, 0.0) * CONST;
		d[7] = (aStage.holes()[hole_num].pointD() + aStage.holes()[hole_num].pointA()) / 2 + Vec2(0.0, 1.0) * CONST;
		for (int i = 0; i < AppNum; i++){
			Vec2 ptoh = d[i] - aStage.player().pos();
			Vec2 htoi = aStage.items()[item_num].pos() - d[i];
			int turn = (ptoh.length() + htoi.length()) / Parameter::PLAYER_SPEED
				+ (Math::Abs(ptoh.rotSign(htoi)) + Math::Abs(playerDir.rotSign(ptoh))) / Parameter::PLAYER_ANGLE_RATE;
			dist_HoleEdge[i] = pair_ii(turn, i);
		}
		qsort_P(dist_HoleEdge, 0, AppNum - 1);
		for (int i = 0; i < AppNum; i++){
			Vec2 v, V;
			v = d[dist_HoleEdge[i].y] - aStage.player().pos();
			//穴に落ちないか
			bool isOK = true;
				Vec2 s1 = aStage.player().pos();
				Vec2 f1 = aStage.player().pos() + v.getNormalized(1.0) * PredictTurn * Parameter::PLAYER_SPEED;
				if (isCross(s1, f1, aStage.holes()[hole_num].pointA(), aStage.holes()[hole_num].pointB())
					|| isCross(s1, f1, aStage.holes()[hole_num].pointB(), aStage.holes()[hole_num].pointC())
					|| isCross(s1, f1, aStage.holes()[hole_num].pointC(), aStage.holes()[hole_num].pointD())
					|| isCross(s1, f1, aStage.holes()[hole_num].pointD(), aStage.holes()[hole_num].pointA())
					)
					isOK=false;
			if (isOK){
				if (Math::IsEqual(aStage.player().vec().cos(v), 1.0f)){
					return Action(ActionType_Move, Parameter::PLAYER_SPEED);
				}
				const float rotSign = aStage.player().vec().rotSign(v);
				// 回転量を絶対値で制限。
				const float rot = Math::LimitAbs(rotSign, Parameter::PLAYER_ANGLE_RATE);
				return Action(ActionType_Rotate, rot);
			}
		}
		std::cout << hole_num << '\n';
		return Action(ActionType_Rotate, 0);
	}

	const int ubTurn = 5; // 4～6
	Action Answer::GetNextAction(const Stage& aStage)
	{
		// 目標にするアイテムを選定。
		//int item_num = -1;
		//greedyに最も近いitem
		int item_num = FindNearestItemNumGreedy(aStage.items(), aStage);
	//	int item_num = FindNearestItemNumGreedyItem(aStage.items(), aStage);
		//convex hull
		//int item_num = FindeNearestItemNumWithConvexHull(aStage.items(), aStage);
		
	//std::ofstream ofs("output.txt", std::ios::app); ofs << currcnt << ' ' << currmax << ' ' << currnum << '\n';
	
		const Vec2 playerDir = aStage.player().vec();
		const Vec2 playerToItem = aStage.items()[item_num].pos() - aStage.player().pos();
		int Turn = Math::Ceil((aStage.items()[item_num].pos() - aStage.player().pos()).length() / Parameter::PLAYER_SPEED);
		const Circle Item(aStage.items()[item_num].pos(), aStage.items()[item_num].radius());
		const Circle PLAYER(aStage.player().pos(), aStage.player().region().radius());
		//ActiveCollisionはそのまま進んだときにItemを取れるかの判定
		const Vec2 Velocity = Parameter::PLAYER_SPEED * playerDir;
		const bool ActiveCollision = WillHit(Item, PLAYER, aStage.player().pos() + Turn * Velocity);
		bool InBoard = aStage.field().isIn(aStage.player().pos() + Turn * Velocity);
		int HoleNum = -1;
		Turn = Math::Min(Turn, ubTurn);
		//そのまま進んで、狙ったアイテムをgetできるか
		//std::ofstream ofs("output.txt", std::ios::app);
		if (ActiveCollision && InBoard) {
			if (aStage.items().count() <=100){
				int currcnt = 0;
				for (int i = 0; i < aStage.items().count(); i++){
					if (aStage.items()[i].isRemoved()) continue;
					const Vec2 playerDir = aStage.player().vec();
					int Turn = Math::Ceil((aStage.items()[item_num].pos() - aStage.player().pos()).length() / Parameter::PLAYER_SPEED);
					const Circle Item(aStage.items()[i].pos(), aStage.items()[i].radius());
					const Circle PLAYER(aStage.player().pos(), aStage.player().region().radius());
					//ActiveCollisionはそのまま進んだときにItemを取れるかの判定
					const Vec2 Velocity = Parameter::PLAYER_SPEED * playerDir;
					const bool ActiveCollision = WillHit(Item, PLAYER, aStage.player().pos() + Turn * Velocity);
					//bool InBoard = aStage.field().isIn(aStage.player().pos() + Turn * Velocity);
					if (ActiveCollision) currcnt++;
				}
				int currnum = 0; int currmax = currcnt;
				for (int k = -3; k <= 3; k++){
					if (k == 0) continue;
					int cnt = 0;
					float ARG = aStage.player().arg() + k * Parameter::PLAYER_ANGLE_RATE;
					const Vec2 playerDir = aStage.player().vec();
					const Circle Item(aStage.items()[item_num].pos(), aStage.items()[item_num].radius());
					for (int i = 0; i < aStage.items().count(); i++){
						if (aStage.items()[i].isRemoved()) continue;
						const Vec2 playerDirTmp = Vec2(Math::Cos(ARG), Math::Sin(ARG));
						const Vec2 playerToItem = aStage.items()[item_num].pos() - aStage.player().pos();
						int Turn = Math::Ceil((aStage.items()[item_num].pos() - aStage.player().pos()).length() / Parameter::PLAYER_SPEED);
						Turn = Math::Min(Turn, ubTurn);
						const Circle ItemTmp(aStage.items()[i].pos(), aStage.items()[i].radius());
						const Circle PLAYER(aStage.player().pos(), aStage.player().region().radius());
						//ActiveCollisionはそのまま進んだときにItemを取れるかの判定
						const Vec2 VelocityTmp = Parameter::PLAYER_SPEED * playerDirTmp;
						const bool ActiveCollisionTmp = WillHit(ItemTmp, PLAYER, aStage.player().pos() + Turn * VelocityTmp);
						const bool ActiveCollision = WillHit(Item, PLAYER, aStage.player().pos() + Turn * VelocityTmp);
						bool InBoard = aStage.field().isIn(aStage.player().pos() + Turn * VelocityTmp);
						if (ActiveCollisionTmp && ActiveCollision && InBoard){
							int HoleNum = -1;
							for (int k = 0; k < aStage.holes().count(); k++){
								VEC playerDir_d = VEC(playerDirTmp.x, playerDirTmp.y);
								VEC S1 = VEC(aStage.player().pos().x, aStage.player().pos().y);
								VEC F1 = S1 + playerDir_d*Turn*(Parameter::PLAYER_SPEED);
								VEC PointA = VEC(aStage.holes()[k].pointA().x, aStage.holes()[k].pointA().y);
								VEC PointB = VEC(aStage.holes()[k].pointB().x, aStage.holes()[k].pointB().y);
								VEC PointC = VEC(aStage.holes()[k].pointC().x, aStage.holes()[k].pointC().y);
								VEC PointD = VEC(aStage.holes()[k].pointD().x, aStage.holes()[k].pointD().y);
								if (isCross(S1, F1, PointA, PointB) || isCross(S1, F1, PointB, PointC)
									|| isCross(S1, F1, PointC, PointD) || isCross(S1, F1, PointD, PointA)
									){
									HoleNum = k; break;
								}
							}
							if (HoleNum == -1) cnt++;
						}
					}
					if (cnt > currmax) currmax = cnt, currnum = k;
				}
				if (currnum > 0) return Action(ActionType_Rotate, Parameter::PLAYER_ANGLE_RATE);
				else if (currnum < 0) return Action(ActionType_Rotate, -1 * Parameter::PLAYER_ANGLE_RATE);
			}

			//ただし、その過程でholeがあればbreak ここは精度高めて効果あり
			for (int j = 0; j < aStage.holes().count(); j++){
				VEC playerDir_d = VEC(playerDir.x, playerDir.y);
				VEC S1 = VEC(aStage.player().pos().x, aStage.player().pos().y);
				VEC F1 = S1 + playerDir_d*Turn*(Parameter::PLAYER_SPEED);
				VEC PointA = VEC(aStage.holes()[j].pointA().x, aStage.holes()[j].pointA().y);
				VEC PointB = VEC(aStage.holes()[j].pointB().x, aStage.holes()[j].pointB().y);
				VEC PointC = VEC(aStage.holes()[j].pointC().x, aStage.holes()[j].pointC().y);
				VEC PointD = VEC(aStage.holes()[j].pointD().x, aStage.holes()[j].pointD().y);
				if (isCross(S1, F1, PointA, PointB) || isCross(S1, F1, PointB, PointC)
					|| isCross(S1, F1, PointC, PointD) || isCross(S1, F1, PointD, PointA)
					){
					HoleNum = j; break;
				}
			}
			//進む過程でどのholeにも交差しないなら、そのまま進む
			//ofs << HoleNum << '\n';
			if (HoleNum == -1) return Action(ActionType_Move, Parameter::PLAYER_SPEED);
			//return SearchWayToEdge(aStage, HoleNum, item_num, Turn);
		}
		bool GOTOHOLE = false; if (HoleNum != -1) GOTOHOLE = true;
		HoleNum = -1;
		for (int i = 0; i < aStage.holes().count(); i++){
			Vec2 s1 = aStage.player().pos();
			Vec2 f1 = aStage.player().pos() + playerToItem.getNormalized(1.0) * Turn * Parameter::PLAYER_SPEED;
			if (isCross(s1, f1, aStage.holes()[i].pointA(), aStage.holes()[i].pointB())
				|| isCross(s1, f1, aStage.holes()[i].pointB(), aStage.holes()[i].pointC())
				|| isCross(s1, f1, aStage.holes()[i].pointC(), aStage.holes()[i].pointD())
				|| isCross(s1, f1, aStage.holes()[i].pointD(), aStage.holes()[i].pointA())
				) {
				HoleNum = i; break;
			}
		}
		//ofs << HoleNum << '\n';
		if (HoleNum != -1)	return SearchWayToEdge(aStage, HoleNum, item_num, Turn);
		
		//理論的には以上で、holeの検出は正しく行われるはずだが、もれがある。
		for (int j = 1; j <= Turn; j++){
			for (int i = 0; i < aStage.holes().count(); i++){
				Vec2 s1 = aStage.player().pos();
				Vec2 f1 = aStage.player().pos() + playerToItem.getNormalized(1.0) * j * Parameter::PLAYER_SPEED;
				if (aStage.holes()[i].isIn(f1)) { HoleNum = i; break; }
			}
		}
		//ofs << HoleNum << '\n';
		if (HoleNum != -1)	return SearchWayToEdge(aStage, HoleNum, item_num, Turn);
		
		if (ActiveCollision && InBoard && !GOTOHOLE) {
			//機能してない
			return Action(ActionType_Move, Parameter::PLAYER_SPEED);
		}
		else {
			// 回転方向付きのなす角を求める。
			const float rotSign = playerDir.rotSign(playerToItem);
			// 回転量を絶対値で制限。
			const float rot = Math::LimitAbs(rotSign, Parameter::PLAYER_ANGLE_RATE);
			return Action(ActionType_Rotate, rot);
		}
	}
}