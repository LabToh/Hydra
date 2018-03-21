/**
*@file CPDecomposition.cpp
*@brief CP分解を行うプログラム
*@author Shogo Satoyama
*@date 2017/01/29
*@details
*CP分解を行うプログラム。N次元に対応している。\n
*hombrew で線形代数ライブラリであるeigenをインストールしておく。\n
*brew install eigen \n
*でインストールできる。
*g++ -O3 -std=c++11 CPDecomposition.cpp -I/usr/local/include \n
*でコンパイルする。
*\n
*入力形式は下記のように数値のあるindexのみにする。\n
*value1, value2, ..., valueN\n
*i, ...\n
*j, ...\n
*k, ...\n
*l, ...\n
*:\n
*と記述\n
*\n
*ex 4階テンソルの例\n
*0.1,0.0,0.11,0.2,0.9\n
*0,1,1,2,3\n
*0,0,0,1,2\n
*0,1,1,3,3\n
*0,3,4,5,9\n
*/

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <limits>
#include <chrono>
#include <random>
#include <iterator>
#include <thread>
#include <random>

using namespace std;

/**
*@brief 画面出力を行う
*@details ベクトルの画面出力を行う。宣言時は、OutPut<double> op;のように型を宣言する
*/
template <typename TYPE>
class OutPut{//出力クラス
public:
  /**
  *@brief コンストラクタ
  */
  OutPut(){};
  /**
  *@brief でコンストラクタ
  */
  virtual ~OutPut(){};
  /**
  *@brief ベクトルを出力
  *@param[in] vector 出力するベクトル
  */
  void PrintVector(vector<TYPE>& vector);//ベクトルを出力
};
template <typename TYPE>
void OutPut<TYPE>::PrintVector(vector<TYPE>& vector){
    for(int NUM = 0; NUM<vector.size(); ++NUM){//ベクトルの数だけループ
        cout << vector[NUM] << ' ';
    }
    cout << endl;
}

/**
*@brief 文字列を分割するクラス
*@details 文字列を分割し、int型、float型、double型、string型のベクトルに格納する
*/
class SplitStep{//1行を分割していく
public:
  /**
  *@brief コンストラクタ
  */
  SplitStep(){};
  /**
  *@brief デコンストラクタ
  */
  virtual ~SplitStep(){};
  /**
  *@brief char型の区切り文字で分割し、vector<int>に格納する
  *@param[out] record 分割によって得られた結果を格納する
  *@param[in] streambuffer 分割する文字列
  *@param[in] delimiter 分割する区切りもじ ','とか
  */
  void Split(vector<int> &record,istringstream &streambuffer,char delimiter);//int型の分割
  /**
  *@brief char型の区切り文字で分割し、vector<float>に格納する
  */
  void Split(vector<float> &record,istringstream &streambuffer,char delimiter);//float型の分割
  /**
  *@brief char型の区切り文字で分割し、vector<double>に格納する
  */
  void Split(vector<double> &record,istringstream &streambuffer,char delimiter);//double型の分割
  /**
  *@brief char型の区切り文字で分割し、vector<string>に格納する
  */
  void Split(vector<string> &record,istringstream &streambuffer,char delimiter);//string型の分割
  /**
  *@brief string型の区切り文字で分割し、vector<string>に格納する
  *@param[out] result 分割によって得られた結果を格納する
  *@param[in] s 分割される文字列
  *@param[in] delim 区切る文字列
  */
  void Splits(vector<string> &result,const std::string& s, const std::string& delim){
    result.clear();
    string::size_type pos = 0;
    while(pos != string::npos ){
      string::size_type p = s.find(delim, pos);
      if(p == string::npos){
        result.push_back(s.substr(pos));
        break;
      }else{
        result.push_back(s.substr(pos, p - pos));
      }
      pos = p + delim.size();
    }
  }
  /**
  *@brief string型の区切り文字で分割し、vector<int>に格納する
  */
  void Splits(vector<int> &result,const std::string& s, const std::string& delim){
    result.clear();
    string::size_type pos = 0;
    while(pos != string::npos ){
      string::size_type p = s.find(delim, pos);
      if(p == string::npos){
        result.push_back(stoi(s.substr(pos)));
        break;
      }else{
        result.push_back(stoi(s.substr(pos, p - pos)));
      }
      pos = p + delim.size();
    }
  }
  /**
  *@brief string型の区切り文字で分割し、vector<float>に格納する
  */
  void Splits(vector<float> &result,const std::string& s, const std::string& delim){
    result.clear();
    string::size_type pos = 0;
    while(pos != string::npos ){
      string::size_type p = s.find(delim, pos);
      if(p == string::npos){
        result.push_back(stof(s.substr(pos)));
        break;
      }else{
        result.push_back(stof(s.substr(pos, p - pos)));
      }
      pos = p + delim.size();
    }
  }
  /**
  *@brief string型の区切り文字で分割し、vector<double>に格納する
  */
  void Splits(vector<double> &result,const std::string& s, const std::string& delim){
    result.clear();
    string::size_type pos = 0;
    while(pos != string::npos ){
      string::size_type p = s.find(delim, pos);
      if(p == string::npos){
        result.push_back(stof(s.substr(pos)));
        break;
      }else{
        result.push_back(stof(s.substr(pos, p - pos)));
      }
      pos = p + delim.size();
    }
  }
};

void SplitStep::Split(vector<int> &record,istringstream &streambuffer,char delimiter){
    string token;// １セル分の文字列
    while (getline(streambuffer, token, delimiter)){// １セル分の文字列をリストに追加する
        record.push_back(stoi(token));//int型に変換
    }
}

void SplitStep::Split(vector<float> &record,istringstream &streambuffer,char delimiter){
    string token;// １セル分の文字列
    while (getline(streambuffer, token, delimiter)){// １セル分の文字列をリストに追加する
        record.push_back(stof(token));//float型に変換
    }
}

void SplitStep::Split(vector<double> &record,istringstream &streambuffer,char delimiter){
    string token;// １セル分の文字列
    while (getline(streambuffer, token, delimiter)){// １セル分の文字列をリストに追加する
        record.push_back(stof(token));//float型に変換
    }
}

void SplitStep::Split(vector<string> &record,istringstream &streambuffer,char delimiter){
    string token;// １セル分の文字列
    while (getline(streambuffer, token, delimiter)){// １セル分の文字列をリストに追加する
        record.push_back(token);//string型のまま
    }
}

/**
*@brief 行ごと、列ごとの和を計算するクラス
*@details Matrixクラスで継承する。
*/
template <typename TYPE>
class Summation{//Σを行うクラス
public:
  /**
  *@brief コンストラクタ
  */
  Summation(){};
  /**
  *@brief デコンストラクタ
  */
  virtual ~Summation(){};
  /**
  *@brief 行ごとの和を計算する
  *@param[out] totalVector 行ごとの和を格納する
  *@param[in] matrix 対象となる行列
  */
  void RowSum(vector<TYPE> &totalVector, vector<vector<TYPE> > matrix){//行ごとの和を計算
    int max_size = matrix.size();
    totalVector.resize(max_size);
    for(int ROW = 0; ROW<max_size; ++ROW){//行の数だけループ
      totalVector[ROW] = accumulate(matrix[ROW].begin(),matrix[ROW].end(),0.0);//行ごとの和をベクトルに格納
    }
  }
  /**
  *@brief 列ごとの和を計算する
  *@param[out] totalVector 列ごとの和を格納する
  *@param[in] matrix 対象となる行列
  */
  void ColSum(vector<TYPE> &totalVector, vector<vector<TYPE> > matrix){//列ごとの和を計算
    TYPE value;
    int max_size = matrix[0].size();
    totalVector.resize(max_size);
    for(int COL = 0; COL<max_size; ++COL){//列の数だけ回す
      value = 0;//初期値
      for(int ROW = 0; ROW<matrix.size(); ++ROW){//行の数だけ回す
        value += matrix[ROW][COL];//合計を計算
      }
      totalVector[COL] = value;//ベクトルに格納
    }
  }
};

/**
*@brief 行列クラス
*@details Matrix<int>、Matrix<string>、Matrix<double>、...で宣言する。\n
*@details ファイルからデータを読み込む時などに使用する。\n
*@details 各行は可変である。
*/
template <typename TYPE>
class Matrix:public OutPut<TYPE>,public SplitStep,public Summation<TYPE>{
public:
  vector<vector<TYPE> > data_;//二次元配列型を格納する箱
  /**
  *@brief コンストラクタ
  */
  Matrix(){};
  /**
  *@brief コンストラクタ
  *@param[in] row_size 行のサイズ
  *@param[in] col_size 列のサイズ
  *@details 確保する行列のサイズがわかっている場合に用いる
  */
  Matrix(int row_size,int col_size){
    data_.resize(row_size);
    SetColSize(col_size);
  }
  /**
  *@brief コンストラクタ
  *@param[in] file_name 読み込むファイルの名前
  *@details CSV形式でデータが記述されている場合はMatrix<double> mat("file_name.csv")とすることでファイルの中身を読み、格納する。
  */
  Matrix(string file_name){
    SetMatrix(file_name,',');
  }
  /**
  *@brief コンストラクタ
  *@param[in] file_name 読み込むファイルの名前
  *@details データがテキストファイルなどに記述されていて、スペースなどでで区切られている場合は、Matrix<double> mat("file_name.csv",' ')とすることでファイルの中身を読み、格納する。
  */
  Matrix(string file_name,char delimiter){
    SetMatrix(file_name,delimiter);
  }
  /**
  *@brief でコンストラクタ
  */
  virtual ~Matrix(){};
  /**
  *@brief 行列の中身の画面表示を行う
  */
  void Print();//表示
  typedef typename vector<TYPE>::const_iterator const_iterator;
  /**
  *@brief ファイルが開けるかの確認
  *@details 開けなければプログラムを終了する
  */
  void OpenCheck(const string& file_name);//ファイルが開けるかの確認
  /**
  *@brief ファイルが開けるかの確認
  *@details 開ければ1、開けなければ-1を返す。
  */
  int CheckFile(const string& file_name);
  /**
  *@brief ファイルからデータを読み込む
  *@param[in] file_name ファイル名
  *@param[in] delimiter 区切り文字
  */
  void SetMatrix(const string& file_name,char delimiter);//CSV形式のファイルを読み込める
  /**
  *@brief ファイルからデータを読み込む
  *@param[in] file_name ファイル名
  *@param[in] delimiter 区切り文字
  *@param[in] skip 読み飛ばす行数を指定
  */
  void SetMatrix(const string& file_name,char delimiter,int skip);//CSV形式のファイルを読み込める
  /**
  *@brief 行の和、列の和を計算する
  *@param[out] vector 計算結果を格納するvector
  *@param[in] option 0で行の和、1で列の和
  *@details R言語のapplyとは違うので気をつける
  */
  void apply(vector<TYPE> &vector, int option);//行ごと列ごとの和の計算をoptionで指定できる
  /**
  *@brief 行列が空か確認する
  *@details 空ならtrue,そうじゃないならfalseが返る
  */
  bool CheckEmpty(){return data_[0].empty();}//行列が空か確認
  /**
  *@brief 行を指定したサイズにリサイズ
  *@param[in] col 列数
  *@param[in]　initial_value 特定の要素で初期化
  */
  void SetColSize(int col,TYPE initial_value){//すべての行を指定した列数にリサイズ
    for(int row = 0; row<data_.size(); ++row){
      data_[row].resize(col,initial_value);
    }
  }
  /**
  *@brief　行を指定したサイズにリサイズ
  *@param[in] col　列数
  */
  void SetColSize(int col){//すべての行を指定した列数にリサイズ
    for(int row = 0; row<data_.size(); ++row){
      data_[row].resize(col);
    }
  }
  /**
  *@brief 特定の行を取得する
  *@param[in] row 取得したい行のindex
  */
  vector<TYPE> const &GetVector (int row) const{return data_[row];}//row行目を取得
  /**
  *@brief 行列積を計算
  *@details C = A*Bといった行列積を計算する場合
  *@details C.data_MatrixProduct(A,B)と入力する。
  */
  void MatrixMatrixProduct(Matrix &matrix1,Matrix& matrix2);//行列の積
  /**
  *@brief ファイルに出力する
  *@param[in]　output_file_name 出力するファイル名
  *@param[in]　delimiter 区切り文字
  */
  void OutPutFile(const string &output_file_name,char delimiter);//ファイルに出力
  /**
  *@brief 特定の行にデータを入れる
  *@param[in] row データを入れたい行
  *@param[in] v 代入するデータ
  */
  void SetVector(int row,vector<TYPE> &v){data_[row] = v;}//row行目にベクトルをセット
  /**
  *@brief　指定した列の、指定した値を見つけvectorにindexを格納する
  *@param[out] index_vector indexを格納するvector
  *@param[in] col 着目する列
  *paraam[in] find_key 着目する値
  */
  void FindAndIndex(vector<int> &index_vector,int col,int find_key);//find_keyを持つ位置を見つけベクトルに格納
  /**
  *@brief　データを削除
  */
  void Clear(){data_.clear();}//matrixを削除
  /**
  *@brief 行列の転置を行う
  *@param[in] input_matrix　転置したい行列
  *@details B = tAとしたい場合
  *@details B.t(A)とする
  */
  void t(Matrix<TYPE> &input_matrix);//転置
  /**
  *@brief　ファイルにあるデータを読み込む
  *@param[in] file_name 読み込むファイル名
  *@param[in] delimiter string型の区切り文字
  *@details　"AAABBBCCC"となっている場合、delimiter = "BBB"とすると、"AAA" "CCC"と区切られる。
  */
  void SetMatrix(const string& file_name,string delimiter){//delimiterを変更することで色々な形式に対応
    OpenCheck(file_name);//ファイルが開けるかチェック
    ifstream filestream(file_name);// ファイルを開く
    while (!filestream.eof()){// ファイルを開く １行読み込む
      string buffer;//string型の箱
      getline(filestream,buffer);// ファイルから読み込んだ１行の文字列を区切り文字で分けてリストに追加する
      if(buffer==""){//最終行などの空白行をSkip
        continue;
      }
      vector<TYPE> record; // １行分の文字列のリスト
      Splits(record,buffer,delimiter);//1行を分割してベクトルに格納
      data_.push_back(record);//行列形式に格納
    }
    filestream.close();
  }
  /**
  *@brief 行列の最大値を要素とするvectorを得る
  *@param[out] vec 各要素が、各行の最大値を格納している。
  */
  void RowMaxVector(vector<TYPE> &vec){//行毎の最大値を要素とするベクトルを入手
    for(int i = 0; i<data_.size(); ++i){
      TYPE temp = *max_element(data_[i].begin(),data_[i].end());//行毎の最大値を取得
      vec.push_back(temp);
    }
  }
  /**
  *@brief 行列全体の最大値を返す
  */
  TYPE MaxElement(){//行列全体の最大値を入手
    vector<TYPE> max_vector;
    RowMaxVector(max_vector);
    TYPE max_value = *max_element(max_vector.begin(),max_vector.end());
    return max_value;
  }
  /**
  *@brief　行毎の最小値を求める
  *@param[in] vec 各要素が、各行の最小値に対応する
  */
  void RowMinVector(vector<TYPE> &vec){//行毎の最小値を要素とするベクトルを入手
    for(int i = 0; i<data_.size(); ++i){
      TYPE temp = *min_element(data_[i].begin(),data_[i].end());//行毎の最大値を取得
      vec.push_back(temp);
    }
  }
  /**
  *@brief　行列全体の最小値を求める
  */
  TYPE MinElement(){//行列全体の最小値を入手
    vector<TYPE> max_vector;
    RowMaxVector(max_vector);
    TYPE min_value = *min_element(max_vector.begin(),max_vector.end());
    return min_value;
  }
};

template <typename TYPE>
void Matrix<TYPE>::Print(){//行列を表示
    cout << "Row = " << data_.size() << endl;//行の数を表示
    cout << "Column = " << data_[0].size() << endl << endl;//列の数を表示
    for(int row = 0; row<data_.size(); ++row){//行の数だけループ
        this->PrintVector(data_[row]);//PrintVectorを呼び出し、1行ごとに表示
    }
    cout << endl;
}
template <class TYPE>
void Matrix<TYPE>::OpenCheck(const string& file_name){//ファイルを開けるかチェック
    fstream filestream(file_name);//ファイルを開く
    if (!filestream.is_open()){//開けない場合
        cout << "SetMatrix, error is occured !!" << endl;//エラーを報告
        cout << file_name <<", Open error !!" << endl;//エラーを報告
        exit(-1);//抜け出る
    }
    filestream.close();//ファイルを閉じる//filestreamを返せるとなお良い
}
template <class TYPE>
int Matrix<TYPE>::CheckFile(const string& file_name){//ファイルを開けるかチェック
    fstream filestream(file_name);//ファイルを開く
    if (!filestream.is_open()){//開けない場合
        cout << "SetMatrix, error is occured !!" << endl;//エラーを報告
        cout << file_name <<", Open error !!" << endl;//エラーを報告
        return -1;
    }
    filestream.close();//ファイルを閉じる//filestreamを返せるとなお良い
    return 1;
}

template <class TYPE>
void Matrix<TYPE>::SetMatrix(const string& file_name,char delimiter){//delimiterを変更することで色々な形式に対応
    OpenCheck(file_name);//ファイルが開けるかチェック
    fstream filestream(file_name);// ファイルを開く
    while (!filestream.eof()){// ファイルを開く １行読み込む
        string buffer;//string型の箱
        getline(filestream,buffer);// ファイルから読み込んだ１行の文字列を区切り文字で分けてリストに追加する
        if(buffer==""){//最終行などの空白行をSkip
          continue;
        }
        vector<TYPE> record; // １行分の文字列のリスト
        istringstream streambuffer(buffer);// 文字列ストリーム
        Split(record,streambuffer,delimiter);//1行を分割してベクトルに格納
        data_.push_back(record);//行列形式に格納
    }
    filestream.close();
}

template <class TYPE>
void Matrix<TYPE>::SetMatrix(const string& file_name,char delimiter,int skip){//delimiterを変更することで色々な形式に対応
    OpenCheck(file_name);//ファイルが開けるかチェック
    fstream filestream(file_name);// ファイルを開く
    int read = 0;
    while (!filestream.eof()){// ファイルを開く １行読み込む
        string buffer;//string型の箱
        getline(filestream,buffer);// ファイルから読み込んだ１行の文字列を区切り文字で分けてリストに追加する
        if(buffer==""){//最終行などの空白行をSkip
            continue;
        }
        if(read<=skip){
            ++read;
            continue;
        }
        vector<TYPE> record; // １行分の文字列のリスト
        istringstream streambuffer(buffer);// 文字列ストリーム
        Split(record,streambuffer,delimiter);//1行を分割してベクトルに格納
        data_.push_back(record);//行列形式に格納
    }
    filestream.close();
}

template <class TYPE>
void Matrix<TYPE>::apply(vector<TYPE> &vector, int option){//行方向、列方向の和を求める
    if(option==0){//optionが０なら行方向
        this->RowSum(vector,data_);//行方向の和を求める
    }
    if(option==1){//optionが1なら列方向
        this->ColSum(vector,data_);//列方向の和を求める
    }
}

template <class TYPE>
void Matrix<TYPE>::MatrixMatrixProduct(Matrix<TYPE> &matrix1,Matrix<TYPE>& matrix2){//行列の積を行う
    if(matrix1.data_[0].size()==matrix2.data_.size()){//積を求める行列の列と行が等しいか確認
        data_.resize(matrix1.data_.size());//格納する行列のサイズを決める
        for(int rowNum = 0; rowNum<matrix1.data_.size(); ++rowNum){//行の数だけ回す
            data_[rowNum].resize(matrix2.data_[0].size());//格納する行列の列のサイズを決める
            for(int element = 0; element<matrix2.data_.size(); ++element){//列の数だけ回す
                if((matrix1.data_[rowNum][element])==0){
                    continue;
                }
                for(int colNum = 0; colNum<matrix2.data_[0].size(); ++colNum){//要素ごとにループ
                    data_[rowNum][colNum] += (matrix1.data_[rowNum][element])*(matrix2.data_[element][colNum]);//積の結果を足していく
                }
            }
        }
    }
else{
        cout <<"Row size " << matrix2.data_.size() << endl;//行と列の値が違った時に行を表示する
        cout <<"Col size " << matrix1.data_[0].size() << endl;//行と列の値が違った時に列を表示する
        cout << "In MatrixMatrixProduct, error is occured !!" << endl;//エラーの報告
        cout << "Row is not equal to column !!" << endl;//エラーの報告
        exit(-1);//エラーが起きたら抜け出る
    }
}

template <class TYPE>
void Matrix<TYPE>::OutPutFile(const string &output_file_name,char delimiter){
    ofstream outputFile(output_file_name);//出力するファイルを開く
    for(int row = 0; row<data_.size(); ++row){//行の数だけ回す
        for(int col = 0; col<data_[row].size(); ++col){//列の数だけ回す
            outputFile << data_[row][col];//区切って出力
            if(col != (data_[row].size()-1)){
                outputFile << delimiter;
            }
        }
        outputFile << endl;//改行を入力
    }
    outputFile.close();//ファイルを閉じる
}

template <class TYPE>
void Matrix<TYPE>::FindAndIndex(vector<int> &index_vector,int col_num,int find_key){
  //col_numにより指定された列中の、find_keyによって指定された要素を持つ列番号を格納する
  for(int row_num = 0; row_num<data_.size(); row_num++){//行の数だけ回す
    if(data_[row_num][col_num]!=find_key){
      continue;//早期continue
    }
    index_vector.push_back(row_num);//find_keyを持つ場合、行番号を格納
  }
}

template <class TYPE>
void Matrix<TYPE>::t(Matrix<TYPE> &input_matrix){//転置を行う
  int set_row_size = input_matrix.data_[0].size();//列の数を取得
  int set_col_size = input_matrix.data_.size();//行の数を取得
  data_.resize(set_row_size);//行を列の数でresize
  SetColSize(set_col_size);//列を行の数でresize
  for(int row = 0; row<input_matrix.data_.size(); ++row){
    for(int col = 0; col<input_matrix.data_[0].size(); ++col){
      data_[col][row] = input_matrix.data_[row][col];//転置を行うi,j -> j,i
    }
  }
}

class RandomEngine{//乱数を生成するクラス
protected:
  vector<double> random_vector_ = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};//0.1刻みの離散値にも対応している
public:
  int RandomInt(int min,int max){//min, maxは範囲に含まれる
    random_device rnd;//非決定的な乱数生成器を生成
    mt19937 mt(rnd());//メルセンヌ・ツイスタの32ビット版、引数は初期シード値
    uniform_int_distribution<> random_int(min, max);//int型のみ乱数で出てくる
    return random_int(mt);//乱数を返す
  }
  double RandomDoubleMinMax(double min, double max){//min, maxは範囲に含まれる
    random_device rnd;//非決定的な乱数生成器を生成
    mt19937 mt(rnd());//メルセンヌ・ツイスタの32ビット版、引数は初期シード値
    uniform_real_distribution<> random_double(min, max);//double型の乱数が出てくる
    return random_double(mt);//乱数を返す
  }
  double RandomDouble(){
    int position = RandomInt(0,10);//indexを生成
    return random_vector_[position];//0.0~1.0の値を返す
  }
  double RandomDouble(int max){
    int position = RandomInt(0,max);//maxは範囲に含まれる
    return random_vector_[position];//0.0~Maxまでの値を返す
  }
};

class TensorDecomposition{//テンソル分解クラス//このクラスを継承しCP分解などを実装
protected:
  vector<Matrix<double> > matrices_;//分解後に作成される行列
public:
  TensorDecomposition(){};
  virtual ~TensorDecomposition(){};
  virtual void Decomposition(int iteration_number){};//テンソル分解
  virtual void Output(string file_name,char delimiter){//ファイルに出力
    for(int i = 0; i<matrices_.size(); ++i){
      string out_file_name = file_name+"_"+to_string(i)+".txt";//ファイル名を作成
      //matrices_[i].Print();
      matrices_[i].OutPutFile(out_file_name,delimiter);//出力
    }
  }
};

class NonNegativeTensorFactorization:public TensorDecomposition{
protected:
  int component_number_;
  vector<int> dimention_vector_;
  string input_file_name_;
  SplitStep spliter_;
  vector<Matrix<double> > one_matrices_;
  vector<double> lambda_vector_;
  double t_error_;//t回目の更新の誤差
  double t_1_error_;//t-1回目の更新の誤差
  double approximation_rate_;//
  double saturation_;
public:
  void SetComponentNumber(int set_component_number){
    component_number_ = set_component_number;
  }
  int GetComponentNumber(){
    return component_number_;
  }
  void SetInputFileName(string set_input_file_name){
    input_file_name_ = set_input_file_name;
  }
  string GetInputFileName(){
    return input_file_name_;
  }
  void SetSaturation(double set_saturation){
    saturation_ = set_saturation;
  }
  double GetSaturation(){
    return saturation_;
  }
  void Initialize(vector<int> dimention_vector){
    matrices_.resize(dimention_vector.size());//リサイズ
    one_matrices_.resize(dimention_vector.size());//リサイズ
    RandomEngine random_engine;//乱数生成機
    for(int m = 0; m<matrices_.size(); ++m){
      matrices_[m].data_.resize(dimention_vector[m]);//行をリサイズ
      one_matrices_[m].data_.resize(dimention_vector[m]);//行をリサイズ
      matrices_[m].SetColSize(component_number_);//列をリサイズ
      one_matrices_[m].SetColSize(component_number_,1.0);//要素を1で初期化
      for(int row = 0; row<dimention_vector[m]; ++row){
        for(int col = 0; col<component_number_; ++col){
          matrices_[m].data_[row][col] = random_engine.RandomDoubleMinMax(0.0,1.0);//初期値として0~1の間で乱数生成
        }
      }
    }
  }
  void InitializeWithNormalize(vector<int> dimention_vector){
    matrices_.resize(dimention_vector.size());//リサイズ
    one_matrices_.resize(dimention_vector.size());//リサイズ
    RandomEngine random_engine;//乱数生成機
    for(int m = 0; m<matrices_.size(); ++m){
      matrices_[m].data_.resize(dimention_vector[m]);//行をリサイズ
      one_matrices_[m].data_.resize(dimention_vector[m]);//行をリサイズ
      matrices_[m].SetColSize(component_number_);//列をリサイズ
      one_matrices_[m].SetColSize(component_number_,1.0);//要素を1で初期化
      for(int row = 0; row<dimention_vector[m]; ++row){
        for(int col = 0; col<component_number_; ++col){
          matrices_[m].data_[row][col] = random_engine.RandomDoubleMinMax(0.0,1.0);//初期値として0~1の間で乱数生成
        }
      }
      Normalize(m);
    }
  }
  double EstimateValue(vector<double> record){//recordには
    double estimated_value = 0.0;//近似によって求められる要素
    int index;
    for(int r = 0; r<component_number_; ++r){//コンポーネントの数だけ回す
      double tmp = 1.0;
      for(int m = 0; m<matrices_.size(); ++m){//分解した行列の数だけ回す
        index = (int)record[(m+1)];//index情報を取得
        tmp *= matrices_[m].data_[index][r];//コンポーネントごとの積を求める
      }
      estimated_value += (tmp);//コンポーネントごとに算出した値をたす
    }
    return estimated_value;
  }
  double EstimateValueWithLambda(vector<double> record){//recordには
    double estimated_value = 0.0;//近似によって求められる要素
    int index;
    for(int r = 0; r<component_number_; ++r){//コンポーネントの数だけ回す
      double tmp = 1.0;
      for(int m = 0; m<matrices_.size(); ++m){//分解した行列の数だけ回す
        index = (int)record[(m+1)];//index情報を取得
        tmp *= matrices_[m].data_[index][r];//コンポーネントごとの積を求める
      }
      estimated_value += (tmp*lambda_vector_[r]);//コンポーネントごとに算出した値をたす
    }
    return estimated_value;
  }
  void Normalize(int matrices_number){
    int row = 0;
    int col = 0;
    vector<double> norm_vector(component_number_,0.0);
    for(col = 0; col<component_number_; ++col){
      for(row = 0; row<matrices_[matrices_number].data_.size(); ++row){
        norm_vector[col] += matrices_[matrices_number].data_[row][col]*matrices_[matrices_number].data_[row][col];
      }
      norm_vector[col] = sqrt(norm_vector[col]);
    }
    for(row = 0; row<matrices_[matrices_number].data_.size(); ++row){
      for(col = 0; col<component_number_; ++col){
        matrices_[matrices_number].data_[row][col] = matrices_[matrices_number].data_[row][col]/norm_vector[col];
      }
    }
    if(matrices_number==0){
      lambda_vector_ = norm_vector;
    }
  }
  void Update(int focus_dimention){
    int index;
    double estimated_value = 0.0;
    double value = 0.0;
    double x;
    fstream filestream(input_file_name_);// ファイルを開く
    Matrix<double> *numerator = new Matrix<double>;//更新式の分子に相当
    Matrix<double> *denominator = new Matrix<double>;//更新式の分母に相当
    numerator->data_ = one_matrices_[(focus_dimention-1)].data_;//要素を1で初期化
    denominator->data_ = one_matrices_[(focus_dimention-1)].data_;//要素を1で初期化
    while (!filestream.eof()){// ファイルを開く １行読み込む
      string buffer;//string型の箱
      getline(filestream,buffer);// ファイルから読み込んだ１行の文字列を区切り文字で分けてリストに追加する
      if(buffer==""){//最終行などの空白行をSkip
        continue;
      }
      vector<double> record; // １行分の文字列のリスト
      istringstream streambuffer(buffer);// 文字列ストリーム
      spliter_.Split(record,streambuffer,',');//1行を分割してベクトルに格納
      value = record[0];//実測値
      estimated_value = EstimateValue(record);
      x = value/estimated_value;//実測値と近似値の比を求める
      //OutPut<double> op;
      //cout << "x = " << x << endl;
      //op.PrintVector(record);
      for(int r = 0; r<component_number_; ++r){//コンポーネントの数だけ回す
        double tmp_value = x;
        double tmp_denominator = 1.0;
        for(int m = 0; m<matrices_.size(); ++m){//分解した行列の数だけ回す
          if(m==(focus_dimention-1)){
            continue;
          }
          index = (int)record[(m+1)];//index情報を取得
          tmp_value = tmp_value*matrices_[m].data_[index][r];//分子部分を計算
          tmp_denominator = tmp_denominator*matrices_[m].data_[index][r];//分母部分を計算
        }
        index = (int)record[focus_dimention];//index情報を取得
        numerator->data_[index][r] = tmp_value;//分子を取得
        denominator->data_[index][r] = tmp_denominator;//分母を取得
      }
    }
    for(int row = 0; row<numerator->data_.size(); ++row){//行の数だけ回す
      for(int col = 0; col<component_number_; ++col){//列の数だけ回す
        matrices_[(focus_dimention-1)].data_[row][col] = matrices_[(focus_dimention-1)].data_[row][col]*(numerator->data_[row][col]/denominator->data_[row][col]);//更新
      }
    }
    filestream.close();//ファイルを閉じる
    delete numerator;
    delete denominator;
  }
  void UpdateWithNormalize(int focus_dimention){
    int index;
    double estimated_value = 0.0;
    double value = 0.0;
    double x;
    fstream filestream(input_file_name_);// ファイルを開く
    Matrix<double> *numerator = new Matrix<double>;//更新式の分子に相当
    Matrix<double> *denominator = new Matrix<double>;//更新式の分母に相当
    numerator->data_ = one_matrices_[(focus_dimention-1)].data_;//要素を1で初期化
    denominator->data_ = one_matrices_[(focus_dimention-1)].data_;//要素を1で初期化
    while (!filestream.eof()){// ファイルを開く １行読み込む
      string buffer;//string型の箱
      getline(filestream,buffer);// ファイルから読み込んだ１行の文字列を区切り文字で分けてリストに追加する
      if(buffer==""){//最終行などの空白行をSkip
        continue;
      }
      vector<double> record; // １行分の文字列のリスト
      istringstream streambuffer(buffer);// 文字列ストリーム
      spliter_.Split(record,streambuffer,',');//1行を分割してベクトルに格納
      value = record[0];//実測値
      estimated_value = EstimateValueWithLambda(record);
      x = value/estimated_value;//実測値と近似値の比を求める
      //OutPut<double> op;
      //cout << "x = " << x << endl;
      //op.PrintVector(record);
      for(int r = 0; r<component_number_; ++r){//コンポーネントの数だけ回す
        double tmp_value = x;
        double tmp_denominator = 1.0;
        for(int m = 0; m<matrices_.size(); ++m){//分解した行列の数だけ回す
          if(m==(focus_dimention-1)){
            continue;
          }
          index = (int)record[(m+1)];//index情報を取得
          tmp_value = tmp_value*matrices_[m].data_[index][r];//分子部分を計算
          tmp_denominator = tmp_denominator*matrices_[m].data_[index][r];//分母部分を計算
        }
        index = (int)record[focus_dimention];//index情報を取得
        numerator->data_[index][r] = tmp_value;//分子を取得
        denominator->data_[index][r] = tmp_denominator;//分母を取得
      }
    }
    for(int row = 0; row<numerator->data_.size(); ++row){//行の数だけ回す
      for(int col = 0; col<component_number_; ++col){//列の数だけ回す
        matrices_[(focus_dimention-1)].data_[row][col] = matrices_[(focus_dimention-1)].data_[row][col]*(numerator->data_[row][col]/denominator->data_[row][col]);//更新
      }
    }
    filestream.close();//ファイルを閉じる
    Normalize((focus_dimention-1));
    delete numerator;
    delete denominator;
  }
  void Update2(int focus_dimention){
    int index;
    double estimated_value = 0.0;
    double value = 0.0;
    double x;
    fstream filestream(input_file_name_);// ファイルを開く
    Matrix<double> *numerator = new Matrix<double>;//更新式の分子に相当
    //Matrix<double> *denominator = new Matrix<double>;//更新式の分母に相当
    numerator->data_ = one_matrices_[(focus_dimention-1)].data_;//要素を1で初期化
    //denominator->data_ = one_matrices_[(focus_dimention-1)].data_;//要素を1で初期化
    while (!filestream.eof()){// ファイルを開く １行読み込む
      string buffer;//string型の箱
      getline(filestream,buffer);// ファイルから読み込んだ１行の文字列を区切り文字で分けてリストに追加する
      if(buffer==""){//最終行などの空白行をSkip
        continue;
      }
      vector<double> record; // １行分の文字列のリスト
      istringstream streambuffer(buffer);// 文字列ストリーム
      spliter_.Split(record,streambuffer,',');//1行を分割してベクトルに格納
      value = record[0];//実測値
      estimated_value = EstimateValue(record);
      x = value/estimated_value;//実測値と近似値の比を求める
      //OutPut<double> op;
      //cout << "x = " << x << endl;
      //op.PrintVector(record);
      for(int r = 0; r<component_number_; ++r){//コンポーネントの数だけ回す
        double tmp_value = x;
        double tmp_denominator = 1.0;
        for(int m = 0; m<matrices_.size(); ++m){//分解した行列の数だけ回す
          if(m==(focus_dimention-1)){
            continue;
          }
          index = (int)record[(m+1)];//index情報を取得
          tmp_value = tmp_value*matrices_[m].data_[index][r];//分子部分を計算
          //tmp_denominator = tmp_denominator*matrices_[m].data_[index][r];//分母部分を計算
        }
        index = (int)record[focus_dimention];//index情報を取得
        numerator->data_[index][r] = tmp_value;//分子を取得
        //denominator->data_[index][r] = tmp_denominator;//分母を取得
      }
    }
    Matrix<double> *denominator = new Matrix<double>;//更新式の分母に相当
    denominator->data_ = one_matrices_[(focus_dimention-1)].data_;//要素を1で初期化
    for(int m2 =0; m2<matrices_.size(); ++m2){
      if((m2+1)==focus_dimention){
        continue;
      }
      vector<double> col_sum;
      matrices_[m2].apply(col_sum,1);
      for(int i = 0; i<denominator->data_.size(); ++i){
        for(int j = 0; j<component_number_; ++j){
          denominator->data_[i][j] = denominator->data_[i][j]/col_sum[j];
        }
      }
    }
    for(int row = 0; row<numerator->data_.size(); ++row){//行の数だけ回す
      for(int col = 0; col<component_number_; ++col){//列の数だけ回す
        matrices_[(focus_dimention-1)].data_[row][col] = matrices_[(focus_dimention-1)].data_[row][col]*(numerator->data_[row][col]/denominator->data_[row][col]);//更新
      }
    }
    filestream.close();//ファイルを閉じる
    delete numerator;
    delete denominator;
  }
  int Check(){
    ComputeError();
    if(abs(t_1_error_ - t_error_)<saturation_){//更新を続けても変化があるか判断
      t_1_error_ = t_error_;//更新
      return 1;
    }else{
      t_1_error_ = t_error_;//更新
      return 0;
    }
  }
  int CheckWithLambda(){
    ComputeErrorWithLambda();
    if(abs(t_1_error_ - t_error_)<saturation_){//更新を続けても変化があるか判断
      t_1_error_ = t_error_;//更新
      return 1;
    }else{
      t_1_error_ = t_error_;//更新
      return 0;
    }
  }
  void ComputeError(){
    int index;
    double estimated_value = 0.0;
    double total_square_estimated_value = 0.0;
    double value = 0.0;
    t_error_ = 0.0;
    OutPut<double> op;
    fstream filestream(input_file_name_);// ファイルを開く
    while (!filestream.eof()){// ファイルを開く １行読み込む
      string buffer;//string型の箱
      getline(filestream,buffer);// ファイルから読み込んだ１行の文字列を区切り文字で分けてリストに追加する
      if(buffer==""){//最終行などの空白行をSkip
        continue;
      }
      vector<double> record; // １行分の文字列のリスト
      istringstream streambuffer(buffer);// 文字列ストリーム
      spliter_.Split(record,streambuffer,',');//1行を分割してベクトルに格納
      value = record[0];
      estimated_value = EstimateValue(record);
      t_error_ += (value-estimated_value)*(value-estimated_value);//二乗誤差
      total_square_estimated_value += estimated_value*estimated_value;
    }
    approximation_rate_ = 1- (t_error_/total_square_estimated_value);
    filestream.close();//ファイルを閉じる
  }
  void ComputeErrorWithLambda(){
    int index;
    double estimated_value = 0.0;
    double total_square_estimated_value = 0.0;
    double value = 0.0;
    t_error_ = 0.0;
    OutPut<double> op;
    fstream filestream(input_file_name_);// ファイルを開く
    while (!filestream.eof()){// ファイルを開く １行読み込む
      string buffer;//string型の箱
      getline(filestream,buffer);// ファイルから読み込んだ１行の文字列を区切り文字で分けてリストに追加する
      if(buffer==""){//最終行などの空白行をSkip
        continue;
      }
      vector<double> record; // １行分の文字列のリスト
      istringstream streambuffer(buffer);// 文字列ストリーム
      spliter_.Split(record,streambuffer,',');//1行を分割してベクトルに格納
      value = record[0];
      estimated_value = EstimateValueWithLambda(record);
      t_error_ += (value-estimated_value)*(value-estimated_value);//二乗誤差
      total_square_estimated_value += estimated_value*estimated_value;
    }
    approximation_rate_ = 1- (t_error_/total_square_estimated_value);
    filestream.close();//ファイルを閉じる
  }
  void RealNormalize(){
    for(int m = 0; m<matrices_.size(); ++m){
      Matrix<double> temp;
      temp.data_ = one_matrices_[m].data_;
      int index;
      double estimated_value = 0.0;
      double value = 0.0;
      OutPut<double> op;
      fstream filestream(input_file_name_);// ファイルを開く
      while (!filestream.eof()){// ファイルを開く １行読み込む
        string buffer;//string型の箱
        getline(filestream,buffer);// ファイルから読み込んだ１行の文字列を区切り文字で分けてリストに追加する
        if(buffer==""){//最終行などの空白行をSkip
          continue;
        }
        vector<double> record; // １行分の文字列のリスト
        istringstream streambuffer(buffer);// 文字列ストリーム
        spliter_.Split(record,streambuffer,',');//1行を分割してベクトルに格納
        value = record[0];
        estimated_value = EstimateValue(record);
        index = (int)record[(m+1)];
        for(int r = 0; r<component_number_; ++r){
          double tijkr = 1.0;
          for(int m2 = 0; m2<matrices_.size(); ++m2){
            int position = (int)record[(m2+1)];
            tijkr = tijkr*matrices_[m2].data_[position][r];
          }
          temp.data_[index][r] += (value/estimated_value)*tijkr;
        }
      }
      filestream.close();//ファイルを閉じる
      temp.OutPutFile("Normalized_"+to_string(m)+".txt",',');
    }
  }
  void Decomposition(int component_number,int iteration_number,vector<int> dimention_vector){
    SetComponentNumber(component_number);
    cout << "Start Non-Negative Tensor Factorization !!" << endl;
    cout << "Start Initialization" << endl;
    Initialize(dimention_vector);
    cout << "Finish Creating Random Matrices." << endl;
    t_1_error_ = 0.0;
    int flag = 0;
    for(int i = 0; i<iteration_number; ++i){
      int focus_dimention = (i%dimention_vector.size()) + 1;//どの行列に着目するかを決定
      Update(focus_dimention);//更新
      /*flag = Check();
      if(flag==1){
        break;
      }*/
      cout << "Iteration: " << i << '\r';
    }
    cout << "Finish Non-Negative Tensor Factorization !!" << endl;
    ComputeError();
    cout << "Squared Error = " << t_error_ << endl;
    cout << "Approximation Rate = " << approximation_rate_ << endl;
  }
  void DecompositionWithLambda(int component_number,int iteration_number,vector<int> dimention_vector){
    SetComponentNumber(component_number);
    cout << "Start Non-Negative Tensor Factorization !!" << endl;
    cout << "Start Initialization" << endl;
    InitializeWithNormalize(dimention_vector);
    cout << "Finish Creating Random Matrices." << endl;
    t_1_error_ = 0.0;
    int flag = 0;
    for(int i = 0; i<iteration_number; ++i){
      int focus_dimention = (i%dimention_vector.size()) + 1;//どの行列に着目するかを決定
      UpdateWithNormalize(focus_dimention);//更新
      /*flag = Check();
      if(flag==1){
        break;
      }*/
      cout << "Iteration: " << i << '\r';
    }
    cout << "Finish Non-Negative Tensor Factorization !!" << endl;
    ComputeErrorWithLambda();
    cout << "Squared Error = " << t_error_ << endl;
    cout << "Approximation Rate = " << approximation_rate_ << endl;
  }
  void Decomposition(int component_number,int iteration_number,vector<int> dimention_vector,double saturation){
    SetComponentNumber(component_number);
    SetSaturation(saturation);
    cout << "Start Non-Negative Tensor Factorization !!" << endl;
    cout << "Start Initialization" << endl;
    Initialize(dimention_vector);
    cout << "Finish Creating Random Matrices." << endl;
    t_1_error_ = 0.0;
    int flag = 0;
    for(int i = 0; i<iteration_number; ++i){
      int focus_dimention = (i%dimention_vector.size()) + 1;//どの行列に着目するかを決定
      Update(focus_dimention);//更新
      flag = Check();
      if(flag==1){
        break;
      }
      cout << "Iteration: " << i << '\r';
    }
    cout << "Finish Non-Negative Tensor Factorization !!" << endl;
    ComputeError();
    cout << "Squared Error = " << t_error_ << endl;
    cout << "Approximation Rate = " << approximation_rate_ << endl;
  }
};

int main(){
  //af::array a = af::constant(1,3,3);
  //af_print(a);
  //af::array b = af::constant(0,3,3);
  //af_print(b);
  NonNegativeTensorFactorization gpu_cp_d;
  gpu_cp_d.SetInputFileName("TestTensorData.txt");
  vector<int> dimention_vector = {10,9,8};
  //gpu_cp_d.SetInputFileName("metagenome_tensor4d.txt");
  //vector<int> dimention_vector = {339,63,31,2};
  //gpu_cp_d.Decomposition(5,1000,dimention_vector,0.0000000001);
  gpu_cp_d.Decomposition(2,10000,dimention_vector);
  gpu_cp_d.RealNormalize();
  gpu_cp_d.Output("output_nonnegative_tensor",',');

  cout << "Tequila Boom Boom !!" << endl;
  return 0;
}
