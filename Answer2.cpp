#include "HPCAnswerInclude.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
// 足し算(float)
float ADD(float a, float b){
	const double EPS = 1e-5;
	if (hpc::Math::Abs(a + b) <= EPS*(hpc::Math::Abs(a) + hpc::Math::Abs(b))) return 0;
	return a + b;
}
// 絶対値(double)
double ABS(double a){
	if (a >= 0) return a; return -a;
}
// 足し算(double)
double ADD(double a, double b){
	const double EPS = 1e-10;
	if (ABS(a + b) <= EPS*(ABS(a) + ABS(b))) return 0;
	return a + b;
}
// STLが使えないので
// pairを準備(float,int),(int,int),(float,float)
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
// このタイプでは実際に演算することがある
struct P{
	float x, y;
	P() {}
	P(float x, float y) : x(x), y(y) {}
	P operator+(P p) { return P(ADD(x, p.x), ADD(y, p.y)); }
	P operator-(P p) { return P(ADD(x, -p.x), ADD(y, -p.y)); }
	P operator*(float d){ return P(x*d, y*d); }
	float dot(P p){ return ADD(x*p.x, y*p.y); }
	float det(P p){ return ADD(x*p.y, -y*p.x); }
};
// double型で同じものを準備
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
	// 新しいステージに入ったときにパラメータなどを
	// 設定し直せるが使っていない
    void Answer::Init(const Stage& aStage)
    {
    }
	// pairに対して、クイックソート
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
	// 円と円の衝突判定(静止しているとき)
	bool IsHit(const Circle& aC0, const Circle& aC1)
	{
		const float squareDist = aC0.pos().squareDist(aC1.pos());
		return squareDist <= (aC0.radius() + aC1.radius()) * (aC0.radius() + aC1.radius());
	}
	// 動いているときの衝突判定(ゲームの問題設定に依存)
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
	// 線分上に点があるのか判定
	bool on_seg(Vec2 s1, Vec2 f1, Vec2 q){
		return ADD((s1 - q).x*(f1 - q).y, -(s1 - q).y*(f1 - q).x) == 0 && ADD((s1 - q).x*(f1 - q).x, (s1 - q).y*(f1 - q).y) <= 0;
	}
	bool on_seg(VEC s1, VEC f1, VEC q){
		return ADD((s1 - q).x*(f1 - q).y, -(s1 - q).y*(f1 - q).x) == 0 && ADD((s1 - q).x*(f1 - q).x, (s1 - q).y*(f1 - q).y) <= 0;
	}
	//線分(s1-f1)と線分(s2-f2)の交点を求める
	Vec2 intersection(Vec2  s1, Vec2 f1, Vec2 s2, Vec2 f2){
		return s1 + (f1 - s1)*((f2 - s2).cross(s2 - s1) / (f2 - s2).cross(f1 - s1));
	}
	VEC intersection(VEC  s1, VEC f1, VEC s2, VEC f2){
		return s1 + (f1 - s1)*((f2 - s2).det(s2 - s1) / (f2 - s2).det(f1 - s1));
	}
	// 二つの線分が交差するのかの判定
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
	// 最も近くにあるアイテムを探す
	const int FindNearestItemNumGreedy(const ItemCollection& items, const Stage& aStage)
	{
		int MIN = 100000; int tmp = 0;
		for (int index = 0; index < items.count(); ++index) {
			// もうすでに取り除かれているアイテムはムシ
		if (items[index].isRemoved()) continue;
			float DIST = (items[index].pos() - aStage.player().pos()).length(); // プレイヤーとアイテムの距離
			Vec2 playerToItem = aStage.items()[index].pos() - aStage.player().pos(); // プレイヤーからアイテムへのベクトル
			Vec2 playerDir = aStage.player().vec(); // プレイヤーの現在の向き
			int MoveTurn = (DIST / Parameter::PLAYER_SPEED); // そのアイテムまで到達するのにかかるターン数
			//衝突判定では(半径r0+r1,中心item.pos())と線分[player.pos(),移動後のplayer.pos()]なので、
			//この円とギリギリ接する点とこの円の中心とplayer.pos()の直角三角形を考えることで、
			//この円との接するのに要するターン数が計算できる
			float l = playerToItem.length(); // プレイヤーとアイテムの距離
			float r = Parameter::PLAYER_RADIUS + items[index].radius();
			//arctan(r/l)
			float arg = Math::ATan2( r, l);

			//ActiveCollisionがあるときは例外
			// 目的としているアイテムの方向に向くまでにかかるターン数
			int RotTurn = ( Math::Abs(aStage.player().vec().rotSign(playerToItem)) - Math::Abs(arg) ) / Parameter::PLAYER_ANGLE_RATE ;
			Vec2 Velocity = Parameter::PLAYER_SPEED * playerDir; // 1ターンで進める度合い(ベクトル)
			// 動いているときの判定をするためにCircle型に代入
			Circle Item(aStage.items()[index].pos(), aStage.items()[index].radius());
			const Circle PLAYER(aStage.player().pos(), aStage.player().region().radius());
			bool ActiveCollision = WillHit(Item, PLAYER, aStage.player().pos() + MoveTurn * Velocity);
			if (ActiveCollision) RotTurn = 0; //そのまま進んでitemを取れる場合は回転する必要がない
			/*
			int cnt = 0;
			if (ActiveCollision){
				RotTurn = 0;
				//積極的に同時に取りに行く
				for (int j = 0; j<items.count(); j++){
					if (items[j].isRemoved()) continue;
					Circle ItemTmp(aStage.items()[j].pos(), aStage.items()[j].radius());
					bool ActiveCollisionTmp = WillHit(ItemTmp, PLAYER, aStage.player().pos() + MoveTurn * Velocity);
					if (ActiveCollisionTmp) cnt++;
				}
			}
			else{
				float ARG = aStage.player().arg();
				if (playerDir.rotSign(playerToItem) >= 0) ARG += Parameter::PLAYER_ANGLE_RATE*RotTurn;
				else ARG -= Parameter::PLAYER_ANGLE_RATE*RotTurn;
				Vec2 playerDirTmp = Vec2(Math::Cos(ARG), Math::Sin(ARG));
				Velocity = Parameter::PLAYER_SPEED * playerDirTmp;
				for (int j = 0; j<items.count(); j++){
					if (items[j].isRemoved()) continue;
					Circle ItemTmp(aStage.items()[j].pos(), aStage.items()[j].radius());
					bool ActiveCollisionTmp = WillHit(ItemTmp, PLAYER, aStage.player().pos() + MoveTurn * Velocity);
					if (ActiveCollisionTmp) cnt++;
				}
				//if (cnt == 0) cnt = 1;
			}
			*/
			int turn = (MoveTurn + RotTurn); // 回転と並進にかかる合計ターン数
			//for (int i = 0; i<cnt; i++) turn *= 0.9;
			if (MIN > turn) MIN = turn, tmp = index; // ターン数最小のものを選択する
		}
		if (MIN != 100000) return tmp; else return -1;
	}
	// 穴に落ちそうな時に穴の周りを通るようにする
	// (穴の端4点を外側に少し拡張してそこを通るようにする)
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
		d[4] = (aStage.holes()[hole_num].pointA() + aStage.holes()[hole_num].pointB()) / 2 + Vec2(-1.0, 0.0) * CONST; // 使用してない
		d[5] = (aStage.holes()[hole_num].pointB() + aStage.holes()[hole_num].pointC()) / 2 + Vec2(0.0, -1.0) * CONST; // 使用してない
		d[6] = (aStage.holes()[hole_num].pointC() + aStage.holes()[hole_num].pointD()) / 2 + Vec2(1.0, 0.0) * CONST; // 使用してない
		d[7] = (aStage.holes()[hole_num].pointD() + aStage.holes()[hole_num].pointA()) / 2 + Vec2(0.0, 1.0) * CONST; // 使用してない
		for (int i = 0; i < AppNum; i++){
			Vec2 ptoh = d[i] - aStage.player().pos(); // プレイヤーから端点へのベクトル
			Vec2 htoi = aStage.items()[item_num].pos() - d[i]; // 端点からアイテムへのベクトル
			// その端点を経由してアイテムを取りに行く場合に要するターン数(回転に要するターン数も含む)
			int turn = (ptoh.length() + htoi.length()) / Parameter::PLAYER_SPEED 
				+ (Math::Abs(ptoh.rotSign(htoi)) + Math::Abs(playerDir.rotSign(ptoh))) / Parameter::PLAYER_ANGLE_RATE;
			dist_HoleEdge[i] = pair_ii(turn, i); // pair(ターン数,端点の番号)
		}
		qsort_P(dist_HoleEdge, 0, AppNum - 1); // クイックソート
		// ターン数が少ない方から見ていく
		for (int i = 0; i < AppNum; i++){
			Vec2 v, V;
			v = d[dist_HoleEdge[i].y] - aStage.player().pos();
			//穴に落ちないか
			bool isOK = true; // trueなら落ちない
			for (int j = 1; j <= PredictTurn; j++){
				V = aStage.player().pos() + v.getNormalized(1.0f)*j*(Parameter::PLAYER_SPEED); // jターン後にプレイヤーのいる場所
				if (aStage.holes()[hole_num].isIn(V)) { isOK = false; break; } // 穴に落ちるならフラグを折る
			}
			if (isOK){
				// プレイヤーの向きと進みたい方向が一致している場合
				if (Math::IsEqual(aStage.player().vec().cos(v), 1.0f)){
					return Action(ActionType_Move, Parameter::PLAYER_SPEED);
				}
				// 回転量を計算
				const float rotSign = aStage.player().vec().rotSign(v);
				// 回転量を絶対値で制限。
				const float rot = Math::LimitAbs(rotSign, Parameter::PLAYER_ANGLE_RATE);
				return Action(ActionType_Rotate, rot);
			}
		}
		return Action(ActionType_Rotate, 0);
	}

	const int ubTurn = 5; // 4～6
    Action Answer::GetNextAction(const Stage& aStage)
    {
        // 目標にするアイテムを選定。
		// とりあえず最も近いアイテムを選ぶ
		int item_num = FindNearestItemNumGreedy(aStage.items(), aStage);
	//std::ofstream ofs("output.txt", std::ios::app); ofs << item_num << '\n';
		const Vec2 playerDir = aStage.player().vec(); // プレイヤーの向き
		const Vec2 playerToItem = aStage.items()[item_num].pos() - aStage.player().pos(); // プレイヤーからアイテムへのベクトル
		// 最も近いアイテムに到達するのに要するターン数
		int Turn = Math::Ceil((aStage.items()[item_num].pos() - aStage.player().pos()).length() / Parameter::PLAYER_SPEED);
		const Circle Item(aStage.items()[item_num].pos(), aStage.items()[item_num].radius());
		const Circle PLAYER(aStage.player().pos(), aStage.player().region().radius());
		//ActiveCollisionはそのまま進んだときにItemを取れるかの判定
		const Vec2 Velocity = Parameter::PLAYER_SPEED * playerDir;
		const bool ActiveCollision = WillHit(Item, PLAYER, aStage.player().pos() + Turn * Velocity);
		bool InBoard = aStage.field().isIn(aStage.player().pos() + Turn * Velocity); // そのまま進んだときにフィールドの中にいるかの判定
		int HoleNum = -1;
		Turn = Math::Min(Turn, ubTurn);
		// そのまま進んで、アイテムをとれて、かつフィールドの中にいる場合
		if (ActiveCollision && InBoard) {
			
			//holeがあればbreak ここは精度高めて効果あり(float->doubleに変更して)
			for (int j = 0; j < aStage.holes().count(); j++){
				VEC playerDir_d = VEC(playerDir.x, playerDir.y);
				VEC S1 = VEC(aStage.player().pos().x, aStage.player().pos().y);
				VEC F1 = S1 + playerDir_d*Turn*(Parameter::PLAYER_SPEED);
				// holeの端点4つ
				VEC PointA = VEC(aStage.holes()[j].pointA().x, aStage.holes()[j].pointA().y);
				VEC PointB = VEC(aStage.holes()[j].pointB().x, aStage.holes()[j].pointB().y);
				VEC PointC = VEC(aStage.holes()[j].pointC().x, aStage.holes()[j].pointC().y);
				VEC PointD = VEC(aStage.holes()[j].pointD().x, aStage.holes()[j].pointD().y);
				// 穴に落ちる場合、どの穴に落ちるかをHoleNumに入れる
				if (isCross(S1, F1, PointA, PointB) || isCross(S1, F1, PointB, PointC)
					|| isCross(S1, F1, PointC, PointD) || isCross(S1, F1, PointD, PointA)
					){
					HoleNum = j; break;
				}
			}
			// 穴に落ちる場合には、それをよけるように行動する
			if (HoleNum != -1) return SearchWayToEdge(aStage, HoleNum, item_num, Turn);
			/*
			if (aStage.items().count() >= 80){
				//現在の方向でitemをいくつ取れるか
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
				//方向を変えてメリットがあれば変えたい
				int currnum = 0; int currmax = currcnt;
				for (int k = -5; k <= 5; k++){
					if (k == 0) continue;
					int cnt = 0;
					float ARG = aStage.player().arg() + k * Parameter::PLAYER_ANGLE_RATE;
					int Turn = Math::Ceil((aStage.items()[item_num].pos() - aStage.player().pos()).length() / Parameter::PLAYER_SPEED);
					const Vec2 playerDir = aStage.player().vec();
					const Vec2 playerDirTmp = Vec2(Math::Cos(ARG), Math::Sin(ARG));
					const Vec2 VelocityTmp = Parameter::PLAYER_SPEED * playerDirTmp;
					const Circle Item(aStage.items()[item_num].pos(), aStage.items()[item_num].radius());
					const bool ActiveCollision = WillHit(Item, PLAYER, aStage.player().pos() + Turn * VelocityTmp);
					if (!ActiveCollision) continue;
					for (int i = 0; i < aStage.items().count(); i++){
						if (aStage.items()[i].isRemoved()) continue;
						const Circle ItemTmp(aStage.items()[i].pos(), aStage.items()[i].radius());
						//ActiveCollisionはそのまま進んだときにItemを取れるかの判定
						const bool ActiveCollisionTmp = WillHit(ItemTmp, PLAYER, aStage.player().pos() + Turn * VelocityTmp);
						bool InBoard = aStage.field().isIn(aStage.player().pos() + Turn * VelocityTmp);
						if (ActiveCollisionTmp && InBoard){
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
			*/
			// 穴に落ちない場合は、そのまま進む
			return Action(ActionType_Move, Parameter::PLAYER_SPEED);
		}
		// 穴に落ちる場合をtrueとする
		bool GOTOHOLE = false; if (HoleNum != -1) GOTOHOLE = true;
		HoleNum = -1;
		// 上では穴の中に落ちる場合を見たが、ここでは
		// 辺(線分)上にある場合を検出する
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
		// 穴に落ちる場合は、よける
		if (HoleNum != -1)	return SearchWayToEdge(aStage, HoleNum, item_num, Turn);
		// そのまま進めばアイテムをとれて、フィールドの中にいて、穴に落ちない場合は、そのまま進む
		if (ActiveCollision && InBoard && !GOTOHOLE) {
			return Action(ActionType_Move, Parameter::PLAYER_SPEED);
		}
		else {
			// 特に何もなければアイテムの方向を向く
            // 回転方向付きのなす角を求める。
            const float rotSign = playerDir.rotSign(playerToItem);
            // 回転量を絶対値で制限。
            const float rot = Math::LimitAbs(rotSign, Parameter::PLAYER_ANGLE_RATE);
            return Action(ActionType_Rotate, rot);
        }
    }
}

//----------------------------------------------------------
// EOF
