#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include <omp.h>
#include <stdexcept>
#include <ranges>
#include <iomanip>
#include <optional>
#include <utility>
#include <algorithm>


class Matrix
{
private:
    std::vector<std::vector<double>> m_mat;
    size_t m_rows;
    size_t m_cols;

public:
    Matrix(int rows,int cols) : m_rows{rows},m_cols{cols}, m_mat{rows,std::vector<double>(cols,0)}  {}
    
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

    Matrix operator* (const Matrix& mat)
    {
        assert(m_cols = mat.m_rows);
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
        return mat;
    }

    std::pair<int,Matrix> gaussElimination() const
    {
        int currentRow = 0;
        int times = 0;
        Matrix copy = *this;
        
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
    
    double determiant()
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

    bool isInvertible()
    {
        return (*this).determiant() != 0;
    }

    Matrix& inverse(Matrix& inverseMat)
    {
        Matrix copy = *this;
        int rows = copy.m_cols;
        Matrix aug{rows,2*rows};

        for (int i=0;i<rows;++i)
        {
            for (int j=0;j<rows;++j)
            {
                aug.m_mat[i][j] = copy.m_mat[i][j];
            }
            aug.m_mat[i][i+rows] = 0.0;
        }

        aug = (aug.gaussElimination()).second;

        for (int k=0;k<rows;++k)
        {
            for (int i = 0;i<rows;++i)
            {
                inverseMat.m_mat[k][i] = aug.m_mat[k][i+rows];
            }
        }
        return inverseMat;
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

    void transpose(Matrix& B,int blocksize)
    {
        for (int i=0;i<m_rows;i+=blocksize)
            for (int j=0;j<m_cols;j+=blocksize)
                for (int ii=i;ii<((i+blocksize > m_rows) ? m_rows : (i+blocksize));++ii)
                    for (int jj=j;jj<((j+blocksize)>m_cols) ? m_cols : (j+blocksize);++jj)
                        B.m_mat[jj][ii] = m_mat[ii][ii];
    }

    double getValue(int r,int c) // 得到mat(r,c)
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

