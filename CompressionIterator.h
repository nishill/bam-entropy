// CompressionIterator.h
// The header interface for the CompressionIterator class

#ifndef __COMPRESSIONITERATOR_H_
#define __COMPRESSIONITERATOR_H_

#include <cstdio> 
#include <string>
#include <stdexcept>
#include <utility>
#include <iostream>
#include <list>
using namespace std;

#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"
using namespace BamTools;

typedef pair<char, char> BaseQualPairFunc;   
typedef map<char, size_t> quality_map;
typedef map<char, quality_map> base_map;
        
class BA_Reader {
    
    private:
        string m_fileName;
        typedef list<base_map> PT;
        typedef base_map BM;
        PT m_pt;
        int32_t m_first_pos;
        void m_insert_base( base_map &, const BamAlignment&, int32_t);    
       
    public:
        
        struct BaseQualInfo {
            char m_base;
            char m_quality;
            size_t m_occurences;
            BaseQualInfo(char base, char quality, size_t occurences ):
                m_base( base ), m_quality( quality ), 
                m_occurences ( occurences ){}
        
        };
        
        class CompressionIterator {

            public:
                CompressionIterator(BA_Reader&);
                CompressionIterator();
                bool operator== ( const CompressionIterator& ) const;
                bool operator!= ( const CompressionIterator& ) const;
                CompressionIterator& operator++ ( int );
                BaseQualInfo operator* ();
            private:
                quality_map::iterator m_qm_iter;
                base_map::iterator m_bm_iter;
                PT::iterator m_pt_iter;
                BA_Reader* m_pba_reader;    
                void m_next();
        };        

        BA_Reader ( const string & fileName );
        void print_tree(); 
        CompressionIterator begin();
        CompressionIterator end(); 
        
        typedef BaseQualPairFunc (*bqpf_summaryFunc)(const base_map& );
        bqpf_summaryFunc m_bqpf_summaryFunc;
        void summarizeBases(bqpf_summaryFunc); 
    
        
        class ListIterator {

            private:
                PT::iterator m_pt_iter;
                bqpf_summaryFunc m_pSummaryFunc;
            public:
                ListIterator( BA_Reader&, bqpf_summaryFunc );
                ListIterator(BA_Reader&);
                bool operator==(const ListIterator&) const;
                bool operator!=(const ListIterator&) const;
                ListIterator& operator++( int );
                pair<char,char> operator*();
        };

        ListIterator lbegin(bqpf_summaryFunc );
        ListIterator lend();    

};

#endif
