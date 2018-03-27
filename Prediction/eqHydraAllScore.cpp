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
      value = 0;
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
    t_row_sum = af::transpose(af::pow(temp_row_sum,lambda_));//λ乗
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

int main(){
  chrono::system_clock::time_point  start, end;
  start = chrono::system_clock::now();
  Matrix<int> first_bipartite("Pathway_Gene.txt");
  Matrix<int> second_bipartite("DiseaseRelatedGene_Pathway.txt");
  Matrix<int> third_bipartite("Disease_DiseaseRelatedGene.txt");
  eqHydra eqHydra(1,0.51944,0.707177);
  eqHydra.first_.bipartite_graph_.TransformMatrixToArrayMatrix(first_bipartite);
  eqHydra.second_.bipartite_graph_.TransformMatrixToArrayMatrix(second_bipartite);
  eqHydra.third_.bipartite_graph_.TransformMatrixToArrayMatrix(third_bipartite);
  eqHydra.first_.bipartite_graph_.EdgeWeighted();
  eqHydra.second_.bipartite_graph_.EdgeWeighted();
  eqHydra.third_.bipartite_graph_.EdgeWeighted();
  cout << "Start eqHydra Prediction." << endl;
  eqHydra.IntegratedAllScore();
  cout << "Processing ..." << endl;
  int N = eqHydra.score_matrix_.array_matrix_.dims(0);
  int M = eqHydra.score_matrix_.array_matrix_.dims(1);
  Matrix<float> result_matrix;
  for(int i = 0; i<N; ++i){
    vector<float> vec(M);
    eqHydra.score_matrix_.array_matrix_.row(i).host(&vec[0]);
    result_matrix.data_.push_back(vec);
  }
  cout << "Finish eqHydra Prediction." << endl;
  result_matrix.OutPutFile("eqHydraPredictionResults.txt",',');
  end = chrono::system_clock::now();
  double elapsed = chrono::duration_cast<chrono::seconds>(end-start).count();
  cout << "Time = " << elapsed << endl;
  cout << "success" << endl;
  return 0;
}
