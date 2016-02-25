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
typedef list<base_map> PT;

// BamAlignmen_Reader class         
class BA_Reader {
    
    private:
        string m_fileName;
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
       
        struct BaseQualCode {
            char m_base;
            char m_quality;
            // coding
            //bool m_inExonRegion;
            BaseQualCode(char base, char quality):
                m_base(base), m_quality(quality){} 
                //m_inExonRegion(inExonRegion){} 
                
        };

        // inner CompressionIterator class 
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
        
        // function pointer utilities 
        //typedef BaseQualPairFunc (*bqpf_summaryFunc)(const base_map& );
        typedef BaseQualCode (*bqpf_summaryFunc)(const base_map&);
        bqpf_summaryFunc m_bqpf_summaryFunc;
        void summarizeBases(bqpf_summaryFunc); 

        char det_mag(char);

        // inner ListIterator class
        class ListIterator {

            private:
                PT::iterator m_pt_iter;
                PT::iterator m_pt_iter_end;
                bqpf_summaryFunc m_pSummaryFunc;
                unsigned m_baseLookaheadBuf[2];    
                unsigned m_qualLookaheadBuf;
                bool m_inExonRegion;
                // start codon Methionine 
                // ATG
                static const unsigned m_Methionine =
                ('A' << 16) | ('T' << 8 ) | ('G' << 0 );
                // stop codon Ochre
                // TAA
                static const unsigned m_Ochre = 
                ('T' << 16) | ('A' << 8 ) | ('A' << 0 );
                // stop codon Amber
                // TAG
                static const unsigned m_Amber =
                ('T' << 16) | ('A' << 8 ) | ('G' << 0 );
                // stop codon Opal
                // TGA
                static const unsigned m_Opal =
                ('T' << 16) | ('G' << 8 ) | ('A' << 0 );
                void m_next();                 
                BaseQualCode m_current();
 
            public:
                ListIterator( BA_Reader&, bqpf_summaryFunc);
                ListIterator(BA_Reader&);
                bool operator==(const ListIterator&) const;
                bool operator!=(const ListIterator&) const;
                ListIterator operator++( int );
                ListIterator& operator++();
                BaseQualCode operator*();
                void m_advance();
                       
        };
        ListIterator lbegin(bqpf_summaryFunc );
        ListIterator lend();    

        
};

#endif

