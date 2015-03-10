HAL_procon2013
==============
HAL研究所プログラミングコンテスト2013の際に使用していたコードです。
HAL研究所プログラミングコンテスト2013の詳細はこちらのurlを参考にしてください。

http://www.hallab.co.jp/progcon/2013/

基本的には、TSP（巡回セールスマン問題)と計算幾何学の問題です。
パラメータの設定上１ターンあたりに回転できる量が少ないので、回転してから進むより
できる限り直線的に進んだ方が効率が良くなります。
当時はプログラミング自体初めて1年ぐらいのときであったので、書きたい処理があまり書けてなかったです。
やったことの概略は以下の通りです。
---
(1) 標準で用意されていたVec2は、floatで定義されていたため、新たにdouble型でVECを定義して使用しました。
これによって誤差のために穴に落ちていた部分がわずかに改善がみられました。

(2) 基本的にはgreedyに近いところを進んでいました。その中で直線的に取れるところはまっすぐに進んでできるだけ
多くのゴミを方向を変えずに取ることを意識して実装しました。このとき、ゴミに対してギリギリ取れるようになることが
多いです。

(3) 進みたい方向に穴がある場合には、穴の端点から長さ√2だけ伸ばした点４つ(例えば、右上の端点から(1,1)だけ離れた点など)の中から、最も近い点を探して通る。

---
