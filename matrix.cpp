#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <ranges>
#include <iomanip>
#include <optional>
#include <utility>
#include <algorithm>
#include <numeric>


class Matrix
{
private:
    std::vector<std::vector<double>> m_mat;
    size_t m_rows;
    size_t m_cols;

public:
    Matrix(int rows,int cols) : m_rows{rows},m_cols{cols}, m_mat{rows,std::vector<double>(cols,0)}  {} // 构造函数，默认全部是0

    // Matrix(Matrix& mat)
    // {
    //     m_mat.resize(mat.m_mat.size());
    //     int i=0;
    //     for (auto& row:m_mat)
    //     {
    //         row.resize(mat.m_mat[0].size());
    //         std::copy(row.begin(),row.end(),mat.m_mat[i].begin());
    //         i++;
    //     }
    // }                  // 析构函数
    
    Matrix operator+ (const Matrix& mat) // 加法
    {
        assert(m_cols == mat.m_cols || m_rows == mat.m_rows);
        Matrix sum_mat{m_rows,m_cols};

        for (size_t r = 0;r < m_rows;++r )
        {
            for (size_t c = 0;c <  m_cols;++c)
            {
                sum_mat.m_mat[r][c] = m_mat[r][c] + mat.m_mat[r][c];
            }
            return sum_mat;
        }

    }

    Matrix operator-(const Matrix& mat)
    {
        assert(m_cols == mat.m_cols || m_rows == mat.m_rows);
        Matrix re_mat{m_rows,m_cols};

        for (size_t r = 0;r < m_rows;++r )
        {
            for (size_t c = 0;c <  m_cols;++c)
            {
                re_mat.m_mat[r][c] = m_mat[r][c] - mat.m_mat[r][c];
            }
            return re_mat;
        }

    }

    bool operator==(const Matrix& other) // 比较两个矩阵是否相等
    {
        if(m_cols != other.m_cols || m_rows != other.m_rows)
        {
            return false;
        }

        for (size_t i=0;i<m_rows;++i)
        {
            if (!std::equal(m_mat[i].begin(),m_mat[i].end(),other.m_mat[i].begin()))
            {
                return false;
            }
        }
        return true;
    }

    bool operator!=(const Matrix& other)
    {
        if(m_cols != other.m_cols || m_rows != other.m_rows)
        {
            return true;
        }

        for (size_t i=0;i<m_rows;++i)
        {
            if (!std::equal(m_mat[i].begin(),m_mat[i].end(),other.m_mat[i].begin()))
            {
                return true;
            }
        }
        return false;
    }
 
    Matrix operator* (const Matrix& mat) const // 矩阵乘法
    {
        if (m_cols != mat.m_rows)
        {
            throw std::invalid_argument("INVALID PARAMETERS");
        }

        Matrix mul_mat{m_rows,mat.m_cols};

        for (size_t r = 0;r<m_rows;++r)
        {
            for (size_t c=0;c<m_cols;++c)
            {
                for (int k = 0;k<m_cols;++k)
                {
                    mul_mat.m_mat[r][c] += m_mat[r][k] * mat.m_mat[k][c];
                }
            }
            
        }
        return mul_mat;
    }

    void rowInterchange(size_t i,size_t j)
    {
        if(i == j || i > m_rows ||j > m_rows)
        {
            throw std::invalid_argument("INVALID PARAMETERS");
        }

        std::swap(m_mat[i],m_mat[j]);
    }

    void multiplication_of_row(size_t row,double factor) // 第二类初等变换
    {
        if (row>m_rows)
        {
            throw std::invalid_argument("INVALID PARAMETERS");
        }

        for (int j=0;j<m_cols;++j)
        {
            m_mat[row][j] *= factor;
        }
    }

    void addition_multiplied_row(size_t row1,size_t row2,double factor)  // 第三类初等变换
    {
        if (row1>m_rows || row2 > m_rows)
        {
            throw std::invalid_argument("INVALID PARAMETERS");
        }

        std::transform(m_mat[row1].begin(),m_mat[row1].begin(),m_mat[row2].begin(),
                                                   m_mat[row1].begin(),[factor](double x,double y){return x + factor*y;});
    }
    

    Matrix subMatrix(size_t start_row,size_t end_row,size_t start_col,size_t end_col ) const //子矩阵确定方法是两条横线，两条竖线交叉产生一个矩阵
    {
        Matrix subMatrix(end_row-start_row,end_col-start_col);
        for (int i = start_row;i<end_col+1;++i)
        {
            std::copy(m_mat[i].begin()+start_col,m_mat[i].begin()+end_col,subMatrix.m_mat[i-start_row].begin());
        }
        return subMatrix;
    }

    std::pair<int,Matrix> gaussElimination() const // 对矩阵进行Gauss消元，返回值中的int times 代表的是进行第三类初等变换的次数，也就是交换两行的次数
    {
        int currentRow = 0;
        int times = 0;
        Matrix copy(*this); // 我们不希望影响矩阵本身
        
        for (int c = 0;c < copy.m_cols && currentRow < m_rows;++c) // traverse every column
        {
             int pRow = currentRow;
             double maxElement = std::abs(copy.m_mat[currentRow][c]);

             for (int r = currentRow+1;r<copy.m_rows;++r)
             {
                if (maxElement < std::abs(copy.m_mat[r][c]))
                {
                    maxElement = std::abs(copy.m_mat[r][c]);
                    pRow = r;
                }
             }

             if (maxElement<1e-10)
             {
                continue;
             }
                    
             
             if (pRow != currentRow)
             {
                std::swap(copy.m_mat[currentRow],copy.m_mat[pRow]);
                times = times + 1;
             } 
            
             

             double pivot = copy.m_mat[currentRow][c];

             for (int r = currentRow+1;r<copy.m_rows;++r)
             {
                double factor = copy.m_mat[r][c] / pivot;
                for (int j = c;j<m_cols;++j)
                {
                    copy.m_mat[r][j] = copy.m_mat[r][j]- factor * copy.m_mat[currentRow][j];
                }

             }

            currentRow ++;
        }
        return {times,std::move(copy)};
    }
    
    double determiant() // 计算行列式，这里需要调用Gauss消元的成员函数
    {
        assert(m_cols == m_rows);
        int det = 1;
        Matrix upTriangle = ((*this).gaussElimination()).second;

        for (int i=0;i<m_cols;++i)
        {
            det = det * upTriangle.m_mat[i][i] * std::pow(-1,((*this).gaussElimination()).first);
        }
        
        return det;
        
    }

    bool isInvertible() // 判断一个矩阵是否是可逆的
    {
        return (*this).determiant() != 0;
    }

    Matrix& inverse(Matrix& inverseMat) const // 
    {
        if(m_cols != m_rows) // 检查是否为方阵
        {
            throw std::invalid_argument("Only square matrix has inverse matrix");
        }

        size_t n = m_rows;
        const double epsilon = 1e-9;

        Matrix augmented_mat(n,2*n); // 构造增广矩阵
        for (size_t i=0;i<n;++i)
        {
            for (size_t j=0;j<n;++j)
            {
                augmented_mat.m_mat[i][j] = m_mat[i][j];
            }
            augmented_mat.m_mat[i][i+n] = 1.0;
        }

        // 我上面写好的gauss消元函数并不能满足使用增广矩阵求逆的需求
        Matrix upper = augmented_mat.gaussElimination().second;
        for (size_t i=0;i<n;++i)
        {
            if (std::abs(upper.m_mat[i][i])<epsilon)
            {
                throw std::runtime_error("This matrix is not invertible");
            }
        }

        Matrix reduced = upper;
        for (int col = n-1;col >= 0;--col) // 从最后一列开始 把增广矩阵的右半部分变为单位矩阵
        {
            double pivot = reduced.m_mat[col][col]; // 定位左边矩阵的主对角线上的元素
            for (size_t j=col;j<2*n;++j) // 主循环，一列一列处理，我们的思路是把主对角线上方的每一个元素都变为0
            {
                reduced.m_mat[col][j] /= pivot; // 先把主对角线上的元素都变为1
            }

            for(size_t row=col-1;row >= 0;--row) // 开始处理主对角线上每个元素上面的所有元素，行索引控制
            {
                double factor = reduced.m_mat[row][col]; // 由于该列主对角线元素已经是1
                for (size_t j = col;j<2*n;++j) // row行其它的元素也要处理
                {
                    reduced.m_mat[row][j] -= factor * reduced.m_mat[col][j];
                }
            }

            // 提取右边的矩阵
            inverseMat = Matrix(n,n);
            for (size_t i =0;i<n;++i)
            {
                for (size_t j=0;j<n;++j)
                {
                    inverseMat.m_mat[i][j] = reduced.m_mat[i][j+n];
                }
            }
            return inverseMat;
        }
        


    }
    
    double trace()
    {
        if (m_cols != m_rows)
        {
            throw std::domain_error("Only square matrix has trace");
        }
        double traceValue{0.0};
        for (int i=0;i<m_cols;++i)
        {
            traceValue = traceValue + m_mat[i][i];
        }
        return traceValue;
    }

    int rank()
    {
        Matrix temp = ((*this).gaussElimination()).second;
        int Rank{1};
        for (int i=0;i<m_rows;++i)
        {
            if (std::all_of(temp.m_mat[i].begin(),temp.m_mat[i].end(),[](int x) { return std::abs(x) <1e-10; }))
            { Rank++;}
        }
        
        return Rank;
    }

    void transpose_by_block(Matrix& transposed,int blocksize)  // 比较适合大型矩阵 但是block的取值是要注意的，可以取二者的gcd
    {
        transposed.m_mat.resize(m_cols);
        for (auto& row:transposed.m_mat)
        {
            row.resize(m_rows);
        }

        for (int i=0;i<m_rows;i+=blocksize)
        {
            for (int j=0;j<m_cols;j+=blocksize)
            {
                // 计算当前块的边界（避免越界）
                int i_end = std::min(i + blocksize, static_cast<int>(m_rows));
                int j_end = std::min(j + blocksize,static_cast<int>(m_cols));

                for (int ii=i;i<i_end;++i)
                {
                    for (int jj=j;jj<j_end;++j)
                    {
                        transposed.m_mat[jj][ii] = m_mat[ii][jj];
                    }
                }
            }
        }
    
    }

    void square_transose()  // 方针的转置
    {
        if (m_cols != m_rows)
        {
            throw std::invalid_argument("INVALID PARAMATERS");
        }

        for (int i=0;i<m_rows;++i)
        {
            for (int j=i+1;j<m_cols;++j)
            {
                std::swap(m_mat[i][j],m_mat[j][i]);
            }
        }        
    }

    double sum() const
    {
        double total = 0.0;
        for (auto& row:m_mat)
        {
            total += std::accumulate(row.begin(),row.end(),0.0);
        }
        return total;
    }

    double product() const
    {
        double total = 1.0;
        for(auto& row:m_mat)
        {
            total *= std::accumulate(row.begin(),row.end(),1.0,[](double acc,double ele){return acc * ele;});
        }

        return total;
    }

    double getValue(int r,int c) const // 得到mat(r,c)
    {
        return m_mat[r][c];
    }
    
    void setValue(int r,int c,double element)
    {
        m_mat[r][c] = element;
    }

    int getCol()
    {
        return m_cols;
    }

    int getRow()
    {
        return m_rows;
    }
    
    void buckFill(const std::vector<std::vector<double>>& data)
    {
        if (data.size() != m_rows && data[0].size() != m_cols)
        {
            throw std::invalid_argument("dismatched dimensions,sorry🥹");
        }

        m_mat = data;
    }

    void identical_fill(double value)
    {
        for (auto& row:m_mat)
        {
            std::fill(row.begin(),row.end(),value);
        }
    }

    void erase() // 清空矩阵
    {
        (*this).identical_fill(0.0);
    }

    double& operator()(int i,int j)
    {
        if(i<1 || i>m_rows || j<1 || j>m_cols)
            throw std::invalid_argument("dismatched dimensions,sorry🥹");

        return m_mat[i-1][j-1];
    }

    friend std::istream& operator>> (std::istream& in, Matrix& mat);
    friend std::ostream& operator<< (std::ostream& out,Matrix& mat);
};

std::istream& operator>> (std::istream& in, Matrix& mat)
    {
        for (int i = 0; i < mat.m_rows; ++i) {
            std::cout << "please input the element of row" << " " << i << '\n';
            for (int j = 0; j < mat.m_cols; ++j) {
                std::cout << "element" << " " << j << ':';
                in >> mat.m_mat[i][j];
            }
        }
        return in;
    }

std::ostream& operator<< (std::ostream& out, Matrix& mat)
{
    std::cout << "the matrix you input is:" << "\n";
    for (int i=0;i<mat.m_rows;++i)
    {
        std::cout << std::left;
        for (int j=0;j<mat.m_cols;++j)
        {
            out << std::setw(10) << mat.m_mat[i][j];
        }
        std::cout << '\n';
    }
    return out;
}

