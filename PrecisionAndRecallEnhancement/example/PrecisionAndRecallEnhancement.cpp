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

class eqHydraRecall{
protected:
  int length_;
public:
  eqHydraRecall(){};
  eqHydraRecall(int length){
    SetLength(length);
  }
  virtual ~eqHydraRecall(){};
  void SetLength(int set_length){
    length_ = set_length;
  }
  void CreateAnswerVector(vector<string> &vec,vector<string> &answer_vec,int answer_position, int recommend_position){
    for(int i = answer_position +1; i<recommend_position; ++i){
      answer_vec.push_back(vec[i]);
    }
  }
  double PointRecall(vector<string> &vec){
    int answer_position = distance(vec.begin(),find(vec.begin(),vec.end(),"Answer:"));
    int recommend_position = distance(vec.begin(),find(vec.begin(),vec.end(),"Recommend:"));
    int score_position = distance(vec.begin(),find(vec.begin(),vec.end(),"Score:"));
    int ignore_position = distance(vec.begin(), find(vec.begin(),vec.end(),"IgnoreIndex:"));
    int end_position = vec.size()-1;
    int number_of_index = score_position-recommend_position-1;
    if(length_>number_of_index){
      cout << "length is larger than number of index " << endl;
      exit(-1);
    }
    vector<string> answer_vec;
    CreateAnswerVector(vec,answer_vec,answer_position,recommend_position);
    double counter = 0;
    for(int i = 0; i<length_; ++i){
      auto itr = find(answer_vec.begin(),answer_vec.end(),vec[recommend_position+1+i]);
      if(itr!=answer_vec.end()){
        ++counter;
      }
    }
    return counter/answer_vec.size();
  }
  double ComputeRecall(Matrix<string> &mat){
    double mean_recall = 0;
    double value = 0;
    for(int i = 0; i<mat.data_.size(); ++i){
      value = PointRecall(mat.data_[i]);
      if(value<0){
        cout << "i = " << i << endl;
        cout << "Something Wrong Recall < 0 " << endl;
        exit(-1);
      }else{
        mean_recall += value;
      }
    }
    mean_recall = mean_recall/mat.data_.size();
    return mean_recall;
  }
  double RecallEnhancement(Matrix<string> &mat,int number_of_diseases){
    double mean_recall = ComputeRecall(mat);
    double recall_enhancement = mean_recall*number_of_diseases/length_;
    return recall_enhancement;
  }
};

class eqHydraPrecision{
protected:
  int length_;
  int links_counter_=0;
public:
  eqHydraPrecision(){};
  eqHydraPrecision(int length){
    SetLength(length);
  }
  void SetLength(int set_length){
    length_ = set_length;
  }
  void CreateAnswerVector(vector<string> &vec,vector<string> &answer_vec,int answer_position, int recommend_position){
    for(int i = answer_position +1; i<recommend_position; ++i){
      answer_vec.push_back(vec[i]);
    }
  }
  double PointPrecision(vector<string> &vec){
    int answer_position = distance(vec.begin(),find(vec.begin(),vec.end(),"Answer:"));
    int recommend_position = distance(vec.begin(),find(vec.begin(),vec.end(),"Recommend:"));
    int score_position = distance(vec.begin(),find(vec.begin(),vec.end(),"Score:"));
    int ignore_position = distance(vec.begin(), find(vec.begin(),vec.end(),"IgnoreIndex:"));
    int end_position = vec.size()-1;
    int number_of_index = score_position-recommend_position-1;
    if(length_>number_of_index){
      cout << "length is larger than number of index " << endl;
      exit(-1);
    }
    vector<string> answer_vec;
    CreateAnswerVector(vec,answer_vec,answer_position,recommend_position);
    links_counter_ += answer_vec.size();
    double counter = 0;
    for(int i = 0; i<length_; ++i){
      auto itr = find(answer_vec.begin(),answer_vec.end(),vec[recommend_position+1+i]);
      if(itr!=answer_vec.end()){
        ++counter;
      }
    }
    return counter/length_;
  }
  double ComputePrecision(Matrix<string> &mat){
    double mean_precision = 0;
    double value = 0;
    links_counter_ = 0;
    for(int i = 0; i<mat.data_.size(); ++i){
      value = PointPrecision(mat.data_[i]);
      if(value<0){
        cout << "i = " << i << endl;
        cout << "Something Wrong Recall < 0 " << endl;
        exit(-1);
      }else{
        mean_precision += value;
      }
    }
    mean_precision = mean_precision/mat.data_.size();
    return mean_precision;
  }
  double PrecisionEnhancement(Matrix<string> &mat,int number_of_diseases,int number_of_genes){
    double mean_precision = ComputePrecision(mat);
    double precision_enhancement = (mean_precision*number_of_diseases*number_of_genes)/links_counter_;
    return precision_enhancement;
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
      double value = -1;
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
  double ComputeAUC(){
    if(totalTP==0&&totalFP==0){
      for(int i = 0; i<score_label_vector.size(); ++i){
        totalTP += score_label_vector[i].second;
      }
      totalFP = score_label_vector.size() -totalTP;
    }
    sort(score_label_vector.begin(),score_label_vector.end());//ソート
    reverse(score_label_vector.begin(),score_label_vector.end());//ソート
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
    return area;
  }
};

int main(int argc, char **argv){
  Matrix<string> mat(argv[1]);
  int length = atoi(argv[2]);
  eqHydraRecall recaller(length);
  eqHydraPrecision precisioner(length);
  ROCScore eqHydra_roc;
  eqHydra_roc.result_matrix_.data_ = mat.data_;
  eqHydra_roc.PrepareROCData();
  int number_of_diseases = 3;
  int number_of_disease_related_genes = 4;
  cout << "Precison Enhancement = " << precisioner.PrecisionEnhancement(mat,number_of_diseases, number_of_disease_related_genes) << endl;
  cout << "Recall Enhancement = " << recaller.RecallEnhancement(mat,number_of_diseases) << endl;
  cout << "AUC = " << eqHydra_roc.ComputeAUC() << endl;
  return 0;
}
