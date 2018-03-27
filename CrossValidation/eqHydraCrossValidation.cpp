#include <arrayfire.h>
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

#define SAMPLES 20
#define ITERATION 100
#include "/path/to/bayesopt/utils/param_loader.hpp"
#include "/path/to/bayesopt/include/bayesopt/bayesopt.hpp"
#include "/path/to/bayesopt/src/bayesoptdisc.cpp"

using namespace std;

template <typename TYPE>
class OutPut{
public:
  OutPut(){};
  virtual ~OutPut(){};
  void PrintVector(vector<TYPE>& vector);
};
template <typename TYPE>
void OutPut<TYPE>::PrintVector(vector<TYPE>& vector){
    for(int NUM = 0; NUM<vector.size(); ++NUM){
        cout << vector[NUM] << ' ';
    }
    cout << endl;
}
class SplitStep{
public:
  SplitStep(){};
  virtual ~SplitStep(){};
  void Split(vector<int> &record,istringstream &streambuffer,char delimiter);
  void Split(vector<float> &record,istringstream &streambuffer,char delimiter);
  void Split(vector<double> &record,istringstream &streambuffer,char delimiter);
  void Split(vector<string> &record,istringstream &streambuffer,char delimiter);
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
    string token;
    while (getline(streambuffer, token, delimiter)){
        record.push_back(stoi(token));
    }
}

void SplitStep::Split(vector<float> &record,istringstream &streambuffer,char delimiter){
    string token;
    while (getline(streambuffer, token, delimiter)){
        record.push_back(stof(token));
    }
}

void SplitStep::Split(vector<double> &record,istringstream &streambuffer,char delimiter){
    string token;
    while (getline(streambuffer, token, delimiter)){
        record.push_back(stof(token));
    }
}

void SplitStep::Split(vector<string> &record,istringstream &streambuffer,char delimiter){
    string token;
    while (getline(streambuffer, token, delimiter)){
        record.push_back(token);
    }
}
template <typename TYPE>
class Summation{
public:
  Summation(){};
  virtual ~Summation(){};
  void RowSum(vector<TYPE> &totalVector, vector<vector<TYPE> > matrix){
    int max_size = matrix.size();
    totalVector.resize(max_size);
    for(int ROW = 0; ROW<max_size; ++ROW){
      totalVector[ROW] = accumulate(matrix[ROW].begin(),matrix[ROW].end(),0.0);
    }
  }
  void ColSum(vector<TYPE> &totalVector, vector<vector<TYPE> > matrix){
    TYPE value;
    int max_size = matrix[0].size();
    totalVector.resize(max_size);
    for(int COL = 0; COL<max_size; ++COL){
      value = 0;//初期値
      for(int ROW = 0; ROW<matrix.size(); ++ROW){
        value += matrix[ROW][COL];
      }
      totalVector[COL] = value;
    }
  }
};

template <typename TYPE>
class Matrix:public OutPut<TYPE>,public SplitStep,public Summation<TYPE>{
public:
  vector<vector<TYPE> > data_;
  Matrix(){};
  Matrix(int row_size,int col_size){
    data_.resize(row_size);
    SetColSize(col_size);
  }
  Matrix(string file_name){
    SetMatrix(file_name,',');
  }
  Matrix(string file_name,char delimiter){
    SetMatrix(file_name,delimiter);
  }
  virtual ~Matrix(){};
  void Print();
  typedef typename vector<TYPE>::const_iterator const_iterator;
  void OpenCheck(const string& file_name);
  int CheckFile(const string& file_name);
  void SetMatrix(const string& file_name,char delimiter);
  void SetMatrix(const string& file_name,char delimiter,int skip);
  void apply(vector<TYPE> &vector, int option);
  bool CheckEmpty(){return data_[0].empty();}
  void SetColSize(int col,TYPE initial_value){
    for(int row = 0; row<data_.size(); ++row){
      data_[row].resize(col,initial_value);
    }
  }
  void SetColSize(int col){
    for(int row = 0; row<data_.size(); ++row){
      data_[row].resize(col);
    }
  }
  vector<TYPE> const &GetVector (int row) const{return data_[row];}
  void MatrixMatrixProduct(Matrix &matrix1,Matrix& matrix2);
  void OutPutFile(const string &output_file_name,char delimiter);
  void SetVector(int row,vector<TYPE> &v){data_[row] = v;}
  void FindAndIndex(vector<int> &index_vector,int col,int find_key);
  void Clear(){data_.clear();}
  void t(Matrix<TYPE> &input_matrix);
  void SetMatrix(const string& file_name,string delimiter){
    OpenCheck(file_name);
    ifstream filestream(file_name);
    while (!filestream.eof()){
      string buffer;
      getline(filestream,buffer);
      if(buffer==""){
        continue;
      }
      vector<TYPE> record;
      Splits(record,buffer,delimiter);
      data_.push_back(record);
    }
    filestream.close();
  }
  void RowMaxVector(vector<TYPE> &vec){
    for(int i = 0; i<data_.size(); ++i){
      TYPE temp = *max_element(data_[i].begin(),data_[i].end());
      vec.push_back(temp);
    }
  }
  TYPE MaxElement(){
    vector<TYPE> max_vector;
    RowMaxVector(max_vector);
    TYPE max_value = *max_element(max_vector.begin(),max_vector.end());
    return max_value;
  }
  void RowMinVector(vector<TYPE> &vec){
    for(int i = 0; i<data_.size(); ++i){
      TYPE temp = *min_element(data_[i].begin(),data_[i].end());
      vec.push_back(temp);
    }
  }
  TYPE MinElement(){
    vector<TYPE> max_vector;
    RowMaxVector(max_vector);
    TYPE min_value = *min_element(max_vector.begin(),max_vector.end());
    return min_value;
  }
};

template <typename TYPE>
void Matrix<TYPE>::Print(){
    cout << "Row = " << data_.size() << endl;
    cout << "Column = " << data_[0].size() << endl << endl;
    for(int row = 0; row<data_.size(); ++row){
        this->PrintVector(data_[row]);
    }
    cout << endl;
}
template <class TYPE>
void Matrix<TYPE>::OpenCheck(const string& file_name){
    fstream filestream(file_name);
    if (!filestream.is_open()){
        cout << "SetMatrix, error is occured !!" << endl;
        cout << file_name <<", Open error !!" << endl;
        exit(-1);
    }
    filestream.close();
}
template <class TYPE>
int Matrix<TYPE>::CheckFile(const string& file_name){
    fstream filestream(file_name);
    if (!filestream.is_open()){
        cout << "SetMatrix, error is occured !!" << endl;
        cout << file_name <<", Open error !!" << endl;
        return -1;
    }
    filestream.close();
    return 1;
}

template <class TYPE>
void Matrix<TYPE>::SetMatrix(const string& file_name,char delimiter){
    OpenCheck(file_name);
    fstream filestream(file_name);
    while (!filestream.eof()){
        string buffer;
        getline(filestream,buffer);
        if(buffer==""){
          continue;
        }
        vector<TYPE> record;
        istringstream streambuffer(buffer);
        Split(record,streambuffer,delimiter);
        data_.push_back(record);
    }
    filestream.close();
}

template <class TYPE>
void Matrix<TYPE>::SetMatrix(const string& file_name,char delimiter,int skip){
    OpenCheck(file_name);
    fstream filestream(file_name);
    int read = 0;
    while (!filestream.eof()){
        string buffer;
        getline(filestream,buffer);
        if(buffer==""){
            continue;
        }
        if(read<=skip){
            ++read;
            continue;
        }
        vector<TYPE> record;
        istringstream streambuffer(buffer);
        Split(record,streambuffer,delimiter);
        data_.push_back(record);
    }
    filestream.close();
}

template <class TYPE>
void Matrix<TYPE>::apply(vector<TYPE> &vector, int option){
    if(option==0){
        this->RowSum(vector,data_);
    }
    if(option==1){
        this->ColSum(vector,data_);
    }
}

template <class TYPE>
void Matrix<TYPE>::MatrixMatrixProduct(Matrix<TYPE> &matrix1,Matrix<TYPE>& matrix2){
    if(matrix1.data_[0].size()==matrix2.data_.size()){
        data_.resize(matrix1.data_.size());
        for(int rowNum = 0; rowNum<matrix1.data_.size(); ++rowNum){
            data_[rowNum].resize(matrix2.data_[0].size());
            for(int element = 0; element<matrix2.data_.size(); ++element){
                if((matrix1.data_[rowNum][element])==0){
                    continue;
                }
                for(int colNum = 0; colNum<matrix2.data_[0].size(); ++colNum){
                    data_[rowNum][colNum] += (matrix1.data_[rowNum][element])*(matrix2.data_[element][colNum]);
                }
            }
        }
    }
else{
        cout <<"Row size " << matrix2.data_.size() << endl;
        cout <<"Col size " << matrix1.data_[0].size() << endl;
        cout << "In MatrixMatrixProduct, error is occured !!" << endl;
        cout << "Row is not equal to column !!" << endl;
        exit(-1);
    }
}

template <class TYPE>
void Matrix<TYPE>::OutPutFile(const string &output_file_name,char delimiter){
    ofstream outputFile(output_file_name);
    for(int row = 0; row<data_.size(); ++row){
        for(int col = 0; col<data_[row].size(); ++col){
            outputFile << data_[row][col];
            if(col != (data_[row].size()-1)){
                outputFile << delimiter;
            }
        }
        outputFile << endl;
    }
    outputFile.close();
}

template <class TYPE>
void Matrix<TYPE>::FindAndIndex(vector<int> &index_vector,int col_num,int find_key){
  for(int row_num = 0; row_num<data_.size(); row_num++){
    if(data_[row_num][col_num]!=find_key){
      continue;
    }
    index_vector.push_back(row_num);
  }
}

template <class TYPE>
void Matrix<TYPE>::t(Matrix<TYPE> &input_matrix){
  int set_row_size = input_matrix.data_[0].size();
  int set_col_size = input_matrix.data_.size();
  data_.resize(set_row_size);
  SetColSize(set_col_size);
  for(int row = 0; row<input_matrix.data_.size(); ++row){
    for(int col = 0; col<input_matrix.data_[0].size(); ++col){
      data_[col][row] = input_matrix.data_[row][col];
    }
  }
}

class ROC :public SplitStep{
  string result_file_name_;
  vector<int> recommend_position_vector_;
  vector<int> number_vector_;
  int length_;
  vector<pair<double,int> > score_label_vector;
  double totalTP = 0;
  double totalFP = 0;
public:
  Matrix<string> result_matrix_;
  void SetLength(int set_length){this->length_ = set_length;}
  void SetFileName(string set_file_name){this->result_file_name_ = set_file_name;}
  ROC(){};
  ROC(string file_name_parameter,int length_parameter){
    SetFileName(file_name_parameter);
    SetLength(length_parameter);
    result_matrix_.SetMatrix(result_file_name_,',');
  }
  void CreatingROCData(){
    for(int row = 0; row<result_matrix_.data_.size(); ++row){
      vector<int> index,ignore;
      vector<double> score_vector;
      int answer_position = distance(result_matrix_.data_[row].begin(),find(result_matrix_.data_[row].begin(),result_matrix_.data_[row].end(),"Answer:"));
      int recommend_position = distance(result_matrix_.data_[row].begin(),find(result_matrix_.data_[row].begin(),result_matrix_.data_[row].end(),"Recommend:"));
      int score_position = distance(result_matrix_.data_[row].begin(),find(result_matrix_.data_[row].begin(),result_matrix_.data_[row].end(),"Score:"));
      int ignore_position = distance(result_matrix_.data_[row].begin(), find(result_matrix_.data_[row].begin(),result_matrix_.data_[row].end(),"IgnoreIndex:"));
      int end_position = result_matrix_.data_[row].size()-1;
      if(end_position == -1){
        break;
      }
      for(int i = recommend_position+1; i<score_position; ++i){
        index.push_back(stoi(result_matrix_.data_[row][i]));
      }
      for(int score_i = score_position+1; score_i<ignore_position; ++score_i){
        score_vector.push_back(stof(result_matrix_.data_[row][score_i]));
      }
      for(int j = ignore_position+1; j<end_position; ++j){
        ignore.push_back(stoi(result_matrix_.data_[row][j]));
      }
      int answer_index = stoi(result_matrix_.data_[row][(answer_position+1)]);
      for(int p = 0; p<index.size(); ++p){
        auto itr = find(ignore.begin(),ignore.end(),index[p]);
        if(itr==ignore.end()){
          if(index[p]==answer_index){
            score_label_vector.push_back(make_pair(score_vector[p],1));
            ++totalTP;
          }else{
            score_label_vector.push_back(make_pair(score_vector[p],0));
            ++totalFP;
          }
        }
      }
    }
  }
  void SetROCData(){
    Matrix<string> temp_matrix;
    cout << "result_file_name_ = " << result_file_name_ << endl;
    temp_matrix.SetMatrix(result_file_name_,',');
    cout <<"matrix size = " <<temp_matrix.data_[0].size() << endl;
    for(int i = 0; i<temp_matrix.data_.size(); ++i){
      score_label_vector.push_back(make_pair(stof(temp_matrix.data_[i][0]),stoi(temp_matrix.data_[i][1])));
    }
  }
  double ComputeAUC(string out_file_name){
    if(totalTP==0&&totalFP==0){
      for(int i = 0; i<score_label_vector.size(); ++i){
        totalTP += score_label_vector[i].second;
      }
      totalFP = score_label_vector.size() -totalTP;
    }
    sort(score_label_vector.begin(),score_label_vector.end());
    reverse(score_label_vector.begin(),score_label_vector.end());
    vector<double> x_pos,y_pos;
    x_pos.push_back(0.0);
    y_pos.push_back(0.0);
    double temp_x = 0.0;
    double temp_y = 0.0;
    for(int j = 0; j<score_label_vector.size(); ++j){
      temp_x += score_label_vector[j].second==0;
      temp_y += score_label_vector[j].second==1;
      x_pos.push_back(temp_x/totalFP);
      y_pos.push_back(temp_y/totalTP);
    }
    double area = 0;
    for(int k = 1; k<x_pos.size(); ++k){
      area += (y_pos[(k-1)]+y_pos[k])*(x_pos[k]-x_pos[(k-1)])/2;
    }
    ofstream ofs(out_file_name);
    for(int l = 1; l<x_pos.size(); ++l){
      ofs << x_pos[l];
      if(l!=(x_pos.size()-1)){
        ofs << ',';
      }
    }
    ofs << endl;
    for(int o = 1; o<y_pos.size(); ++o){
      ofs << y_pos[o];
      if(o!=(y_pos.size()-1)){
        ofs << ',';
      }
    }
    ofs.close();
    return area;
  }
};

class ROCScore{
public:
  Matrix<string> result_matrix_;
  string result_file_name_;
  vector<pair<double,int> > score_label_vector;
  double totalTP = 0;
  double totalFP = 0;
public:
  ROCScore(){};
  virtual ~ROCScore(){};
  void SetFileName(string set_file_name){this->result_file_name_ = set_file_name;}
  ROCScore(string file_name_parameter){
    SetFileName(file_name_parameter);
    result_matrix_.SetMatrix(result_file_name_,',');
  }
  int PositionSearch(vector<string> vec,string key_word){
    int position = distance(vec.begin(),find(vec.begin(),vec.end(),key_word));
    return position;
  }
  void CreateVector(vector<string> input,int pos1,int pos2,vector<int> &output){
    for(int i = pos1+1; i<pos2; ++i){
      output.push_back(stoi(input[i]));
    }
  }
  void CreateVector(vector<string> input,int pos1,int pos2,vector<double> &output){
    for(int i = pos1+1; i<pos2; ++i){
      output.push_back(stof(input[i]));
    }
  }
  void PrepareROCData(){
    for(int row = 0; row<result_matrix_.data_.size(); ++row){
      vector<int> index_vector,ignore_vector,answer_vector;
      vector<double> score_vector;
      int answer_position = PositionSearch(result_matrix_.data_[row],"Answer:");
      int recommend_position = PositionSearch(result_matrix_.data_[row],"Recommend:");
      int score_position = PositionSearch(result_matrix_.data_[row],"Score:");
      int ignore_position = PositionSearch(result_matrix_.data_[row],"IgnoreIndex:");
      int end_position = result_matrix_.data_[row].size()-1;
      if(end_position == -1){
        break;
      }
      CreateVector(result_matrix_.data_[row],answer_position,recommend_position,answer_vector);
      CreateVector(result_matrix_.data_[row],recommend_position,score_position,index_vector);
      CreateVector(result_matrix_.data_[row],score_position,ignore_position,score_vector);
      CreateVector(result_matrix_.data_[row],ignore_position,end_position,ignore_vector);
      double value = -1;//scoreは絶対に負の値を取らないため
      int rank = 0;
      int answer_index = stoi(result_matrix_.data_[row][(answer_position+1)]);
      for(int answer_i = answer_position+1; answer_i<recommend_position; ++answer_i){
        answer_vector.push_back(stoi(result_matrix_.data_[row][answer_i]));
      }
      for(int p = 0; p<index_vector.size(); ++p){
        auto itr = find(answer_vector.begin(),answer_vector.end(),index_vector[p]);
        if(itr!=answer_vector.end()){
          score_label_vector.push_back(make_pair(score_vector[p],1));
          ++totalTP;
        }else{
          score_label_vector.push_back(make_pair(score_vector[p],0));
          ++totalFP;
        }
      }
    }
  }
  void SetROCData(){
    Matrix<string> temp_matrix;
    temp_matrix.SetMatrix(result_file_name_,',');
    for(int i = 0; i<temp_matrix.data_.size(); ++i){
      score_label_vector.push_back(make_pair(stof(temp_matrix.data_[i][0]),stoi(temp_matrix.data_[i][1])));
    }
  }
  double ComputeAUC(string out_file_name){
    if(totalTP==0&&totalFP==0){
      for(int i = 0; i<score_label_vector.size(); ++i){
        totalTP += score_label_vector[i].second;
      }
      totalFP = score_label_vector.size() -totalTP;
    }
    sort(score_label_vector.begin(),score_label_vector.end());
    reverse(score_label_vector.begin(),score_label_vector.end());
    vector<double> x_pos,y_pos;
    double temp_x = 0.0;
    double temp_y = 0.0;
    for(int j = 0; j<score_label_vector.size(); ++j){
      temp_x += score_label_vector[j].second==0;
      temp_y += score_label_vector[j].second==1;
      x_pos.push_back(temp_x/totalFP);
      y_pos.push_back(temp_y/totalTP);
    }
    double area = 0;
    for(int k = 0; k<x_pos.size(); ++k){
      if(k==0){
        area += (y_pos[k])*(x_pos[k])/2;
      }else{
        area += (y_pos[(k-1)]+y_pos[k])*(x_pos[k]-x_pos[(k-1)])/2;
      }
    }
    Matrix<double> mat,t_mat;
    mat.data_.push_back(x_pos);
    mat.data_.push_back(y_pos);
    t_mat.t(mat);
    t_mat.OutPutFile(out_file_name,',');
    return area;
  }
};

class Remover{
public:
  Remover(){};
  virtual ~Remover(){};
  void RemoveEmptyCols(af::array &ar,af::array &result){
    af::array t = anyTrue(ar,0);
    af::array col_ind = where(t);
    result = ar(af::span, col_ind);
  }
  void RemoveEmptyRows(af::array &ar,af::array &result){
    af::array t = anyTrue(ar,1);
    af::array row_ind = where(t);
    result =  ar(row_ind, af::span);
  }
};

class ArrayMatrix:public Remover{
public:
  af::array array_matrix_;
  ArrayMatrix(){};
  ArrayMatrix(int row_size,int col_size){
    array_matrix_ = af::array(row_size,col_size);
  }
  virtual ~ArrayMatrix(){};
  void Initialize(int row_size,int col_size){
    array_matrix_ = af::array(row_size,col_size);
  }
  void TransformMatrixToArrayMatrix(Matrix<int> &bipartite_graph){
    int row_size = bipartite_graph.data_.size();
    int col_size = bipartite_graph.data_[0].size();
    array_matrix_ = af::array(row_size,col_size);
    for(int i = 0; i<row_size; ++i){
      vector<int> vec = bipartite_graph.data_[i];
      af::array temp_array(1,col_size,vec.data());
      array_matrix_.row(i) = temp_array;
    }
  }
  void TransformMatrixToArrayMatrix(Matrix<float> &bipartite_graph){
    int row_size = bipartite_graph.data_.size();
    int col_size = bipartite_graph.data_[0].size();
    array_matrix_ = af::array(row_size,col_size);
    for(int i = 0; i<row_size; ++i){
      vector<float> vec = bipartite_graph.data_[i];
      af::array temp_array(1,col_size,vec.data());
      array_matrix_.row(i) = temp_array;
    }
  }
  void TransformMatrixToArrayMatrix(Matrix<int> &bipartite_graph,af::array &array_mat){
    int row_size = bipartite_graph.data_.size();
    int col_size = bipartite_graph.data_[0].size();
    array_mat = af::array(row_size,col_size);
    for(int i = 0; i<row_size; ++i){
      vector<int> vec = bipartite_graph.data_[i];
      af::array temp_array(1,col_size,vec.data());
      array_mat.row(i) = temp_array;
    }
  }
  void EdgeWeighted(){
    af::array row_sum,col_sum;
    col_sum = af::sum(array_matrix_,0);
    array_matrix_ = array_matrix_/tile(col_sum,array_matrix_.dims(0));
  }
};

class Hybrid{
protected:
  double lambda_;
public:
  ArrayMatrix bipartite_graph_;
  ArrayMatrix score_matrix_;
  ArrayMatrix weight_matrix_;
  Hybrid(){};
  virtual ~Hybrid(){};
  Hybrid(float lambda_parameter){
    SetLambda(lambda_parameter);
  }
  void SetLambda(float set_lambda){
    lambda_ = set_lambda;
  }
  float GetLambda(){
    return lambda_;
  }
  virtual void AllScore(){
    RowComputingWeight();
    score_matrix_.array_matrix_ = af::matmul(weight_matrix_.array_matrix_,bipartite_graph_.array_matrix_);
  }
  void RowComputingWeight(){
    af::array row_sum,temp_row_sum,t_row_sum,col_sum,pow_matrix,col_sum_matrix;
    temp_row_sum = af::sum(bipartite_graph_.array_matrix_,1);
    row_sum = af::pow(temp_row_sum,1-lambda_);
    t_row_sum = af::transpose(af::pow(temp_row_sum,lambda_));
    pow_matrix = af::matmul(row_sum,t_row_sum);
    col_sum = af::sum(bipartite_graph_.array_matrix_,0);
    col_sum_matrix = tile(col_sum,bipartite_graph_.array_matrix_.dims(0));
    weight_matrix_.array_matrix_ = af::matmul(bipartite_graph_.array_matrix_/col_sum_matrix,af::transpose(bipartite_graph_.array_matrix_))/pow_matrix;
  }
  void ColComputingWeight(){
    af::array computing_array_matrix = transpose(bipartite_graph_.array_matrix_);
    af::array row_sum_0,col_sum_0, col_sum_matrix;
    row_sum_0 = af::sum(computing_array_matrix,1);
    af::dim4 row_sum_1_dim(row_sum_0.dims(0),1);
    af::dim4 row_sum_2_dim(1,row_sum_0.dims(0));
    af::array row_sum_1 = af::moddims(af::pow(row_sum_0,1-lambda_),row_sum_1_dim);
    af::array row_sum_2 = af::moddims(af::pow(row_sum_0,lambda_),row_sum_2_dim);
    af::array pow_matrix = af::matmul(row_sum_1,row_sum_2);
    col_sum_0 = af::sum(computing_array_matrix,0);
    col_sum_matrix = tile(col_sum_0,computing_array_matrix.dims(0));
    weight_matrix_.array_matrix_ = af::matmul(computing_array_matrix/col_sum_matrix,bipartite_graph_.array_matrix_)/pow_matrix;
  }
  void RowComputingEdgeWeightedWeight(){
    bipartite_graph_.EdgeWeighted();
    RowComputingWeight();
  }
  void ColComputingEdgeWeightedWeight(){
    bipartite_graph_.EdgeWeighted();
    ColComputingWeight();
  }
};

class CrossValidation{
public:
  Matrix<int> index_matrix_;
  CrossValidation(){};
  virtual ~CrossValidation(){};
  void CreateIndexing(Matrix<int> &bipartite_graph){
    vector<int> row_degree, col_degree;
    bipartite_graph.apply(row_degree,0);
    bipartite_graph.apply(col_degree,1);
    int n = bipartite_graph.data_.size();
    int m = bipartite_graph.data_[0].size();
    vector<int> temp_index;
    temp_index.resize(2);
    for(int i = 0; i<n; ++i){
      if(row_degree[i]==1){
        continue;
      }
      for(int j = 0; j<m; ++j){
        if(bipartite_graph.data_[i][j]==0){
          continue;
        }else{
          if(col_degree[j]==1){
            continue;
          }
          temp_index[0] = i;
          temp_index[1] = j;
          index_matrix_.data_.push_back(temp_index);
        }
      }
    }
  }
  virtual void Process(){};
  virtual void Jackknife(){
    for(int i = 0; i<index_matrix_.data_.size(); ++i){
      Process();
    }
  }
};

class eqHydra{
public:
  Hybrid first_;
  Hybrid second_;
  Hybrid third_;
  ArrayMatrix score_matrix_;
  eqHydra(){};
  eqHydra(float first_lambda,float second_lambda,float third_lambda){
    first_.SetLambda(first_lambda);
    second_.SetLambda(second_lambda);
    third_.SetLambda(third_lambda);
  }
  virtual ~eqHydra(){};
  void IntegratedAllScore(){
    first_.AllScore();
    second_.AllScore();
    third_.AllScore();
    score_matrix_.array_matrix_ = af::matmul(third_.score_matrix_.array_matrix_,af::matmul(second_.score_matrix_.array_matrix_,first_.score_matrix_.array_matrix_));
  }
};

class GeneSelection{
public:
  vector<int> gene_selection_vector_;
  GeneSelection(){};
  virtual ~GeneSelection(){};
  void Selection(Matrix<int> third_phase,Matrix<int> second_phase){
    Matrix<int> temp_third,temp_second;
    for(int gene = 0; gene<third_phase.data_[0].size(); ++gene){
      cout << "gene = " << gene << endl;
      temp_second.data_ = second_phase.data_;
      temp_third.data_ = third_phase.data_;
      for(int i = 0; i<third_phase.data_.size(); ++i){
        if(temp_third.data_[i][gene]==1){
          temp_third.data_[i][gene] = 0;
        }
      }
      vector<int> third_row_sum;
      temp_third.apply(third_row_sum,0);
      auto itr_third = find(third_row_sum.begin(),third_row_sum.end(),0);
      if(itr_third!=third_row_sum.end()){
        continue;
      }
      for(int j = 0; j<second_phase.data_[0].size(); ++j){
        if(temp_second.data_[gene][j]==1){
          temp_second.data_[gene][j] = 0;
        }
      }
      vector<int> second_col_sum;
      temp_second.apply(second_col_sum,1);
      auto itr_second = find(second_col_sum.begin(),second_col_sum.end(),0);
      if(itr_second!=second_col_sum.end()){
        continue;
      }
      if(itr_third==third_row_sum.end()&&itr_second==second_col_sum.end()){
        gene_selection_vector_.push_back(gene);
      }
    }
  }
};

class eqHydraCrossValidation:public Remover, public GeneSelection,public SplitStep{
protected:
  string out_file_name;
public:
  vector<int> index_seq;
  vector<float> score_seq;
  Matrix<int> third_bipartite_graph_;
  Matrix<int> disease_gene_link_matrix_;
  af::array preserve_second_;
  af::array preserve_third_;
  eqHydra eqHydra_;
  eqHydraCrossValidation(){};
  virtual ~eqHydraCrossValidation(){};
  void SetFileName(string set_out_file_name){
    out_file_name = set_out_file_name;
  }
  string GetFileName(){
    return out_file_name;
  }
  void Preserve(){
    preserve_third_ = eqHydra_.third_.bipartite_graph_.array_matrix_;
    preserve_second_ = eqHydra_.second_.bipartite_graph_.array_matrix_;
  }
  void CreateDiseaseLink(){
    disease_gene_link_matrix_.data_.resize(third_bipartite_graph_.data_[0].size());
    for(int disease = 0; disease<third_bipartite_graph_.data_.size(); ++disease){
      for(int disease_gene = 0; disease_gene<third_bipartite_graph_.data_[0].size(); ++disease_gene){
        if(third_bipartite_graph_.data_[disease][disease_gene]==0){
          continue;//コンティニュー
        }else{//リンクがあれば
          disease_gene_link_matrix_.data_[disease_gene].push_back(disease);
        }
      }
    }
  }
  virtual void Process(ROCScore &roc,int i){
    int position = gene_selection_vector_[i];
    eqHydra_.second_.bipartite_graph_.array_matrix_ = preserve_second_;
    eqHydra_.third_.bipartite_graph_.array_matrix_ = preserve_third_;
    int N = eqHydra_.third_.bipartite_graph_.array_matrix_.dims(0);
    int L = eqHydra_.second_.bipartite_graph_.array_matrix_.dims(1);
    gfor(af::array a_i,N){
      eqHydra_.third_.bipartite_graph_.array_matrix_(a_i,position) = 0;
    }
    gfor(af::array a_j,L){
      eqHydra_.second_.bipartite_graph_.array_matrix_(position,a_j) = 0;
    }
    RemoveEmptyCols(eqHydra_.third_.bipartite_graph_.array_matrix_,eqHydra_.third_.bipartite_graph_.array_matrix_);
    RemoveEmptyRows(eqHydra_.second_.bipartite_graph_.array_matrix_,eqHydra_.second_.bipartite_graph_.array_matrix_);
    eqHydra_.second_.bipartite_graph_.EdgeWeighted();
    eqHydra_.second_.AllScore();
    eqHydra_.third_.AllScore();
    eqHydra_.score_matrix_.array_matrix_ = af::matmul(eqHydra_.third_.score_matrix_.array_matrix_,af::matmul(eqHydra_.second_.score_matrix_.array_matrix_,eqHydra_.first_.score_matrix_.array_matrix_.col(position)));
    vector<string> record;
    int flag = StringConcatenation(record,position,eqHydra_.score_matrix_.array_matrix_);
    if(flag==1){
      roc.result_matrix_.data_.push_back(record);
    }
    af::deviceGC();
  };
  virtual double Jackknife(float first_lambda,float second_lambda,float third_lamda){
    eqHydra_.first_.SetLambda(first_lambda);
    eqHydra_.second_.SetLambda(second_lambda);
    eqHydra_.third_.SetLambda(third_lamda);
    eqHydra_.first_.AllScore();
    ROCScore *roc = new ROCScore;
    for(int i = 0; i<gene_selection_vector_.size(); ++i){
      Process(*roc,i);
    }
    string temp_out_file_name = out_file_name+"FirstLambda"+to_string(first_lambda)+"SecondLambda"+to_string(second_lambda)+"ThirdLambda"+to_string(third_lamda)+".txt";
    roc->result_matrix_.OutPutFile(temp_out_file_name,',');
    string roc_file_name = "ROC"+temp_out_file_name;
    roc->PrepareROCData();
    double value = roc->ComputeAUC(roc_file_name);
    delete roc;
    return value;
  }
  int StringConcatenation(vector<string> &record,int j,af::array &score){
    int N = score.dims(0);
    af::array out;
    af::array indice;
    af::sort(out,indice,score,0,false);
    index_seq.resize(N);
    score_seq.resize(N);
    out.host(&score_seq[0]);
    indice.host(&index_seq[0]);
    record.push_back("FocusUser:");
    record.push_back(to_string(j));
    record.push_back("Answer:");
    for(int answer = 0; answer<disease_gene_link_matrix_.data_[j].size(); ++answer){
      record.push_back(to_string(disease_gene_link_matrix_.data_[j][answer]));
    }
    record.push_back("Recommend:");
    for(int element = 0; element<N; ++element){
      record.push_back(to_string(index_seq[element]));
    }
    record.push_back("Score:");
    for(int score = 0; score<N; ++score){
      record.push_back(to_string(score_seq[score]));
    }
    record.push_back("IgnoreIndex:");
    record.push_back("End");
    index_seq.clear();
    score_seq.clear();
    return 1;
  }
};

class eqHydraBayesianOptimization: public bayesopt::ContinuousModel{
 public:
  eqHydraCrossValidation eqhydra_cv;
  void OutFileNameSetting(string set_out_file_name){
    eqhydra_cv.SetFileName(set_out_file_name);
  }
  void BipartiteSetting(Matrix<int> &first_matrix,Matrix<int> &second_matrix,Matrix<int> &third_matrix,Matrix<int> &selection_matrix){
    for(int i = 0; i<selection_matrix.data_[0].size(); ++i){
      eqhydra_cv.gene_selection_vector_.push_back(selection_matrix.data_[0][i]);
    }
    eqhydra_cv.third_bipartite_graph_.data_ = third_matrix.data_;
    eqhydra_cv.CreateDiseaseLink();
    eqhydra_cv.eqHydra_.first_.bipartite_graph_.TransformMatrixToArrayMatrix(first_matrix);
    eqhydra_cv.eqHydra_.second_.bipartite_graph_.TransformMatrixToArrayMatrix(second_matrix);
    eqhydra_cv.eqHydra_.third_.bipartite_graph_.TransformMatrixToArrayMatrix(third_matrix);
    eqhydra_cv.eqHydra_.first_.bipartite_graph_.EdgeWeighted();
    eqhydra_cv.eqHydra_.third_.bipartite_graph_.EdgeWeighted();
    eqhydra_cv.Preserve();
  }
  void BipartiteSetting(Matrix<int> &first_matrix,Matrix<int> &second_matrix,Matrix<int> &third_matrix){
    eqhydra_cv.Selection(third_matrix,second_matrix);
    eqhydra_cv.third_bipartite_graph_.data_ = third_matrix.data_;
    eqhydra_cv.CreateDiseaseLink();
    eqhydra_cv.eqHydra_.first_.bipartite_graph_.TransformMatrixToArrayMatrix(first_matrix);
    eqhydra_cv.eqHydra_.second_.bipartite_graph_.TransformMatrixToArrayMatrix(second_matrix);
    eqhydra_cv.eqHydra_.third_.bipartite_graph_.TransformMatrixToArrayMatrix(third_matrix);
    eqhydra_cv.eqHydra_.first_.bipartite_graph_.EdgeWeighted();
    eqhydra_cv.eqHydra_.third_.bipartite_graph_.EdgeWeighted();
    eqhydra_cv.Preserve();
  }
  eqHydraBayesianOptimization(bayesopt::Parameters param):
  ContinuousModel(3,param) {}
  double evaluateSample( const boost::numeric::ublas::vector<double> &query ){
    double value = eqhydra_cv.Jackknife(query[0],query[1],query[2]);
    return (-1.0)*value;
  };
  bool checkReachability(const vectord &query)
  {return true;};
};


int main(){
  chrono::system_clock::time_point  start, end;
  start = chrono::system_clock::now();
  Matrix<int> first_bipartite("Pathway_Gene.txt");
  Matrix<int> second_bipartite("DiseaseRelatedGene_Pathway.txt");
  Matrix<int> third_bipartite("Disease_DiseaseRelatedGene.txt");
  Matrix<int> selection_matrix("DiseaseGeneIndexForCrossValidation.txt");

  bayesopt::Parameters parameters;
  parameters.l_type = L_EMPIRICAL;
  parameters.n_init_samples = SAMPLES;
  parameters.n_iterations = ITERATION;
  eqHydraBayesianOptimization eqhydrabo(parameters);
  int parameter_number = 3;
  vectord result(parameter_number);
  vectord upper(parameter_number);
  vectord lower(parameter_number);
  for(int element = 0; element<parameter_number; ++element){
    upper[element] = 1.0;
    lower[element] = 0.0;
  }
  eqhydrabo.setBoundingBox(lower,upper);
  eqhydrabo.BipartiteSetting(first_bipartite,second_bipartite,third_bipartite,selection_matrix);
  eqhydrabo.OutFileNameSetting("eqHydra(pathway)");
  eqhydrabo.optimize(result);
  end = chrono::system_clock::now();
  double elapsed = chrono::duration_cast<chrono::seconds>(end-start).count();
  cout << "Time = " << elapsed << endl;
  cout << "success" << endl;
  return 0;
}
