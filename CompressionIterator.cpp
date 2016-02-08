// CompressionIterator.cpp
// A class implementation of CompressionIterator.h

#include <cstdio>
#include <cstring>
#include <cctype>
#include <iomanip>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <stdexcept>
#include "CompressionIterator.h"
using namespace std;

#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"
#include "api/SamSequence.h"
#include "api/SamSequenceDictionary.h"
using namespace BamTools;

void BA_Reader::print_ops() {
    int32_t curr_pos = m_first_pos;
    PT::iterator pos_iter = m_pt.begin();
    CigarOp op = '\0';
    while ( pos_iter != m_pt.end()) {
        cout << "Genome position " << curr_pos  << endl;
        base_map::iterator base_iter = pos_iter->begin();
        while ( base_iter != pos_iter->end()) {
            char curr_base = base_iter->first;
            cout << "\t Current base is " << curr_base <<endl;
            base_iter++;
        }
        pos_iter++;
        curr_pos++;
    }
    
}

// print_tree()
void BA_Reader::print_tree() {
    
    int32_t curr_pos = m_first_pos;
    PT::iterator pos_iter = m_pt.begin();
    while ( pos_iter != m_pt.end()) {
        cout << "Genome position " << curr_pos  << endl;
        base_map::iterator base_iter = pos_iter->begin();
        while ( base_iter != pos_iter->end()) {
            char curr_base = base_iter->first;
            cout << "\t Current base is " << curr_base <<endl;
            quality_map::iterator qual_iter = base_iter->second.begin();
            while (qual_iter != base_iter->second.end()) {
                cout << "\t\t Quality scores are " << 
                qual_iter->first<< " => "  << qual_iter->second << endl;
                qual_iter++;
            }
            base_iter++;
        }
        pos_iter++;
        curr_pos++;
    }

} 

// m_insert_base()
void BA_Reader::m_insert_base( base_map & p_map, 
                    const BamAlignment & ba, int32_t base_pos ){ 
    char base_key = ba.AlignedBases[base_pos];
    char quality_key = ba.Qualities[base_pos];
    // indicates inconsistency within the bam file
    if ( !isgraph(base_key) or !isgraph(quality_key)) return;
    
    typedef pair < base_map :: iterator, bool > BaseResult;
    BaseResult base_result =
            p_map.emplace(base_key, quality_map());
    typedef pair < quality_map :: iterator, bool > QualityResult;
    QualityResult qr = 
            base_result.first->second.emplace(quality_key, 0u);
    qr.first->second++;
}

// BamAlignmentReader() constructor
BA_Reader::BA_Reader ( const string & fileName ) : 
    m_fileName ( fileName ){

    BamReader bar;
    bar.Open( fileName );
    BamAlignment alignment;
    vector<CigarOp> cigar_data = alignment.CigarData;
    PT::iterator pnext_map = m_pt.end();
    PT::iterator pcurr_map = m_pt.end();
    bool success = bar.GetNextAlignment(alignment);
    if (success ) m_first_pos = alignment.Position;
    while ( success ) {
        int32_t curr_pos = alignment.Position;
        
        /*vector<CigarOp>::iterator vco = cigar_data.begin();
        for (; vco != cigar_data.end(); vco++ ) {
            //cout << vco->Length;
            for ( int i = 0; i < vco->Length; i++ ) {
                cout << *i << endl;
            }
        }*/
        BamAlignment nextAlignment;
        success = bar.GetNextAlignment ( nextAlignment );
        int32_t next_pos;
        if ( success ) { 
            next_pos = nextAlignment.Position;
        }
        else {
            next_pos = -1;
        }
        for ( int32_t i = 0; i < alignment.Length; i++ ) {
            if ( pcurr_map == m_pt.end() ) {
                pcurr_map = m_pt.emplace(m_pt.end());
            } 
            m_insert_base( *pcurr_map, alignment, i );
            if ( i + curr_pos == next_pos ) {
                pnext_map = pcurr_map;   
            }
            pcurr_map++;
        }
        pcurr_map = pnext_map;
        pnext_map = m_pt.end();
        alignment = nextAlignment;
    }

} 

// BamAlignmentReader::begin()
BA_Reader::CompressionIterator BA_Reader::begin() {
    return CompressionIterator (*this);
}

// BamAlignmentReader::end()
BA_Reader::CompressionIterator BA_Reader::end(){
    return CompressionIterator(); 
}

// CompressionIterator::operator==()
bool BA_Reader::CompressionIterator::operator==
        ( const CompressionIterator& rhs) const { 
    if ( m_pba_reader != rhs.m_pba_reader ) {
        return false;
    }
    else if (m_pba_reader == NULL ) {
        return true;
    }
    else {
        return m_pt_iter == rhs.m_pt_iter and
               m_pba_reader == rhs.m_pba_reader and
               m_bm_iter == rhs.m_bm_iter;   

    }
}

// CompressionIterator::operator!=()
bool BA_Reader::CompressionIterator::operator!=
        ( const CompressionIterator& rhs) const {
    return !this->operator== (rhs);
}

// BA_Reader::CompressionIterator() constructor
BA_Reader::CompressionIterator::CompressionIterator(BA_Reader& reader):
    m_pt_iter(reader.m_pt.begin()), m_pba_reader(& reader) {
    while ( m_pt_iter != m_pba_reader->m_pt.end() ) {
        m_bm_iter = m_pt_iter->begin ();
        while ( m_bm_iter != m_pt_iter->end() ) {
            m_qm_iter = m_bm_iter->second.begin ();
            if ( m_qm_iter != m_bm_iter->second.end() ) {
                return;
            }
            m_bm_iter++;
        }
        m_pt_iter++;
    }
    m_pba_reader = NULL;
}
// constructor without arguments 
BA_Reader::CompressionIterator::CompressionIterator():
    m_pba_reader(NULL) {}


// BA_Reader::ComressionIterator::m_next()
void BA_Reader::CompressionIterator::m_next(){
    while( m_pt_iter != m_pba_reader->m_pt.end()) {
        while ( m_bm_iter != m_pt_iter->end()) {
            if ( m_qm_iter != m_bm_iter->second.end()) {
                return;
            }
            m_bm_iter++;
            if (m_bm_iter != m_pt_iter->end()) {
                m_qm_iter = m_bm_iter->second.begin();
            }
        }
        m_pt_iter++;
        if (m_pt_iter != m_pba_reader->m_pt.end()) {
            m_bm_iter = m_pt_iter->begin();
        }
    }
    m_pba_reader = NULL;
}

// operator++(int)
BA_Reader::CompressionIterator & BA_Reader::
    CompressionIterator::operator++( int ){
    m_qm_iter++;
    m_next ();
    return *this;
}
 
// operator*()
BA_Reader::BaseQualInfo BA_Reader::CompressionIterator::operator * (){
    return BaseQualInfo ( m_bm_iter->first, m_qm_iter->first,
                            m_qm_iter->second );
}


// summarizeBases()
void BA_Reader::summarizeBases(bqpf_summaryFunc pFunc ) {
    PT::iterator iter = m_pt.begin();
    while ( iter != m_pt.end() ) {
        BaseQualPairFunc bqpf = (*pFunc)(*iter);
        char base = bqpf.first;
        char qual = bqpf.second;
        iter++;
    }
}



//////////////////////////////////////
////    ListIterator Class        ////
//////////////////////////////////////

// ListIterator constructor for begin
BA_Reader::ListIterator::ListIterator(BA_Reader& read, 
                                        bqpf_summaryFunc sumFunc): 
    m_pt_iter(read.m_pt.begin()), m_pSummaryFunc(sumFunc){}

// ListIterator constructor for end
BA_Reader::ListIterator::ListIterator(BA_Reader& read ):
    m_pt_iter(read.m_pt.end()), m_pSummaryFunc(NULL) {}  

// ListIterator lbegin()
BA_Reader::ListIterator BA_Reader::lbegin(bqpf_summaryFunc p_func) {
    return ListIterator(*this, p_func);
}

// ListIterator lend()
BA_Reader::ListIterator BA_Reader::lend(){
    return ListIterator(*this);
}

// ListIterator operator++
BA_Reader::ListIterator& BA_Reader::ListIterator::operator++ ( int ) {
    BA_Reader::ListIterator temp = *this;
    m_pt_iter++;
    return temp; 
}

// ListIterator operator++
BA_Reader::ListIterator& BA_Reader::ListIterator::operator++() {
    m_pt_iter++;
    return *this;
}

// ListIterator operator*
pair<char,char> BA_Reader::ListIterator::operator*() {
    pair <char,char> summary = (*m_pSummaryFunc)(*m_pt_iter);
    assert(isgraph(summary.first));
    assert(isgraph(summary.second));
    return summary;
}

// ListIterator operator==
bool BA_Reader::ListIterator::operator==(const ListIterator& rhs) const{
    return m_pt_iter == rhs.m_pt_iter;
}

// ListIterator operator!=
bool BA_Reader::ListIterator::operator!=(const ListIterator& rhs ) const {
    return !this->operator==(rhs);
}
