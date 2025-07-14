#include <linearList.h>
#include <iostream>

template <typename T>
struct  matrixTerm
{
   int row;
   int col;
   T value;
};

template <typename T>
class spareMatrix
{
private:
    arrayList<matrixTerm<T>> m_terms;
    int m_rows;
    int m_cols;

public:
    int rows() const { return m_rows; }
    int cols() const { return m_cols; }


    void transpose(spareMatrix<T>& b) // 输出参数 实现稀疏矩阵的快速转置
    {
        b.m_cols = m_rows;
        b.m_rows = m_cols;
        b.m_terms.reset(m_terms.size());

        int* colSize = new int[m_cols+1];
        int* rowNext = new int[m_cols+1];

        for (int i=1;i<= m_cols;++i)
        {
            col[i] = 0;
        }
        
        for (typename arrayList<matrixTerm<T>>::iterator it = m_terms.begin();it != m_terms.end();++it)
        {
            ++colSize[it->col];
        }

        rowNext[1] = 0;
        for (int i=2;i<=m_cols;++i)
        {
            rowNext[i] = rowNext[i-1] + colSize[i-1];
        }

        matrixTerm<T> temp;
        for (typename arrayList<matrixTerm<T>>::iterator it = m_terms.begin();it != m_terms.end();++it)
        {
            int j = rowNext[it->col]++;
            temp.col = it->row;
            temp.row = it->col;
            temp.value = it->value;
            b.m_terms.set(j,temp);
        }
    }

    void add(spareMatrix<T>&b,spareMatrix<T>& c) // 稀疏矩阵的相加
    {
        if(m_rows != b.m_rows || m_cols != b.m_cols)
        {
            throw std::invalid_argument();
        }

        c.m_cols = m_cols;
        c.m_rows = m_rows;
        c.m_terms.clear();
        int cSize{0};

        typename arrayList<matrixTerm<T>>::iterator it = m_terms.begin();
        typename arrayList<matrixTerm<T>>::iterator ib = b.m_terms.begin();
        typename arrayList<matrixTerm<T>>::iterator itEnd = m_terms.end();
        typename arrayList<matrixTerm<T>>::iterator ibEnd = b.m_terms.end();

        while (it!= itEnd && ib!= ibEnd)
        {
            int itIndex = (it->row)*m_cols + (it->col);
            int ibIndex = (ib->row)*m_cols + (ib->col);

            if (itIndex<ibIndex)
            {
                c.m_terms.insert(cSize++,*it);
                it++;
            }

            if (itIndex==ibIndex)
            {
                if((it->value) + (ib->value) != 0 )
                {
                    matrixTerm<T> temp;
                    temp.col = it->col;
                    temp.row = it->row;
                    temp.value = (it->value) + (ib->value);
                    c.m_terms.insert(cSize++,temp);
                }

                ++it;
                ++ib;
            }

            else
            {
                c.m_terms.insert(cSize++,*ib);
                ++ib;
            }
        }

        for (;it!=itEnd;++it)
        {
            c.m_terms.insert(cSize++,*it);
        }

        for (;ib!=ibEnd;++ib)
        {
            c.m_terms.insert(cSize++,*ib);
        }
    }
    
    
    
    

};
template <typename T>
std::ostream& operator<<(std::ostream& out,spareMatrix<T>& x) 
{
    out << "row: " << x.rows() << " " << "columns: " << x.cols() << std::endl;
    out << "nonzero terms: " << x.m_terms.size() << std::endl;

    for (typename arrayList<matrixTerm<T> >::iterator it = x.m_terms.begin();it != x.m_terms.end();++it)
    {
        out << "a(" << it->row << ", " << it->col << ") = " << it->value << std::endl;
    }

    return out;
}

template <typename T>
std::istream& operator>> (std::istream& in,spareMatrix<matrixTerm<T>>& x)
{
    int number;
    std::cout << "enter rows,cols and number of terms" << std::endl;
    in >> x.rows() >> x.cols() << number ;
    if (x.rows() < 0 || x.cols() <0 || number <0)
    {
        std::cout << "输入无效，请输入三个positive数！" << std::endl;
        in.clear();  // 清除错误标志
        in.ignore(10000, '\n');  // 忽略缓冲区剩余内容

        std::cout << "enter rows,cols and number of terms" << std::endl;
        in >> x.rows() >> x.cols() << number ;
    }
    
    x.m_terms.reset(number);
    matrixTerm<T> mTerm;

    for(int i=0;i<number;++i)
    {
        std::cout << "enter data here row col value" << i+1 << std::endl;
        in >> mTerm.row >> mTerm.col << mTerm.value;

        x.m_terms.set(i,mTerm);
    }

    return in;
}
