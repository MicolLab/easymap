
/*
 * Copyright 2018, Chanhee Park <parkchanhee@gmail.com> and Daehwan Kim <infphilo@gmail.com>
 *
 * This file is part of HISAT 2.
 *
 * HISAT 2 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HISAT 2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HISAT 2.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "timer.h"
#include "aligner_sw.h"
#include "aligner_result.h"
#include "scoring.h"
#include "sstring.h"

#include "repeat_builder.h"
#include "repeat_kmer.h"
#include "bit_packed_array.h"

unsigned int levenshtein_distance(const std::string& s1, const std::string& s2)
{
    const std::size_t len1 = s1.size(), len2 = s2.size();
    std::vector<unsigned int> col(len2+1), prevCol(len2+1);

    for (unsigned int i = 0; i < prevCol.size(); i++)
        prevCol[i] = i;
    for (unsigned int i = 0; i < len1; i++) {
        col[0] = i+1; 
        for (unsigned int j = 0; j < len2; j++)
            // note that std::min({arg1, arg2, arg3}) works only in C++11,
            // for C++98 use std::min(std::min(arg1, arg2), arg3)
            col[j+1] = std::min( std::min(prevCol[1 + j] + 1, col[j] + 1), prevCol[j] + (s1[i]==s2[j] ? 0 : 1) );
        col.swap(prevCol);
    }
    return prevCol[len2];
}


unsigned int hamming_distance(const string& s1, const string& s2)
{
    assert_eq(s1.length(), s2.length());

    uint32_t cnt = 0;

    for(size_t i = 0; i < s1.length(); i++) {
        if(s1[i] != s2[i]) {
            cnt++;
        }
    }

    return cnt;
}

string reverse(const string& str)
{
    string rev = str;
    size_t str_len = str.length();
    
    for(size_t i = 0; i < str_len; i++) {
        rev[i] = str[str_len - i - 1];
    }
    
    return rev;
}

string reverse_complement(const string& str)
{
    string rev = str;
    size_t str_len = str.length();
    
    for(size_t i = 0; i < str_len; i++) {
        char nt = str[str_len - i - 1];
        rev[i] = asc2dnacomp[nt];
    }
    
    return rev;
}

template<typename TStr>
TIndexOffU getLCP(const TStr& s,
                  CoordHelper& coordHelper,
                  TIndexOffU a,
                  TIndexOffU b,
                  TIndexOffU max_len = 1000)
{
    TIndexOffU a_end = coordHelper.getEnd(a);
    TIndexOffU b_end = coordHelper.getEnd(b);
    
    assert_leq(a_end, s.length());
    assert_leq(b_end, s.length());
    
    TIndexOffU k = 0;
    while((a + k) < a_end && (b + k) < b_end && k < max_len) {
        if(s[a + k] != s[b + k]) {
            break;
        }
        k++;
    }
    
    return k;
}

template<typename TStr>
bool isSameSequenceUpto(const TStr& s,
                              CoordHelper& coordHelper,
                              TIndexOffU a,
                              TIndexOffU b,
                              TIndexOffU upto)
{
    TIndexOffU a_end = coordHelper.getEnd(a);
    TIndexOffU b_end = coordHelper.getEnd(b);
    assert_leq(a_end, s.length());
    assert_leq(b_end, s.length());
    
    if(a + upto > a_end || b + upto > b_end)
        return false;

    for(int k = upto - 1; k >= 0; k--) {
        if(s[a + k] != s[b + k]) {
            return false;
        }
    }
    
    return true;
}

bool isSenseDominant(CoordHelper& coordHelper,
                     const EList<TIndexOffU>& positions,
                     size_t seed_len)
{
    // count the number of seeds on the sense strand
    size_t sense_mer_count = 0;
    for(size_t i = 0; i < positions.size(); i++) {
        if(positions[i] < coordHelper.forward_length())
            sense_mer_count++;
    }
    
    // there exists two sets of seeds that are essentially the same
    // due to sense and antisense strands except palindromic sequences
    // choose one and skip the other one
    assert_leq(sense_mer_count, positions.size());
    size_t antisense_mer_count = positions.size() - sense_mer_count;
    if(sense_mer_count < antisense_mer_count) {
        return false;
    } else if(sense_mer_count == antisense_mer_count) {
        assert_geq(positions.back(), coordHelper.forward_length());
        TIndexOffU sense_pos = coordHelper.length() - positions.back() - seed_len;
        if(positions[0] > sense_pos) return false;
    }
    
    return true;
}

static const Range EMPTY_RANGE = Range(1, 0);

struct RepeatRange {
	RepeatRange() {
        forward = true;
    };
	RepeatRange(Range r, int id) : 
		range(r), rg_id(id), forward(true) {};
    RepeatRange(Range r, int id, bool fw) :
        range(r), rg_id(id), forward(fw) {};

	Range range;
	int rg_id;
    bool forward;
};

Range reverseRange(const Range& r, TIndexOffU size)
{
	size_t len = r.second - r.first;
	Range rc;

	rc.first = size - r.second;
	rc.second = rc.first + len;

	return rc;
}

string reverseComplement(const string& str)
{
	string rc;
	for(TIndexOffU si = str.length(); si > 0; si--) {
		rc.push_back(asc2dnacomp[str[si - 1]]);
	}
	return rc;
}


string toMDZ(const EList<Edit>& edits, const string& read)
{
    StackedAln stk;
    BTDnaString btread;
    BTString buf;

    btread.install(read.c_str(), true);

    stk.init(btread, edits, 0, 0, 0, 0);
    stk.buildMdz();
    stk.writeMdz(&buf, NULL);

    return string(buf.toZBuf());
}

string applyEdits(const string& ref, size_t read_len, const EList<Edit>& edits, const Coord& coord)
{
    string read;
    size_t ref_pos = coord.off();
    size_t read_pos = 0;

    for(size_t i = 0; i < edits.size(); i++) {
        for(; read_pos < edits[i].pos; read_pos++, ref_pos++) {
            read.push_back(ref[ref_pos]);
        }

        if(edits[i].type == EDIT_TYPE_READ_GAP) {
            // delete on read
            ref_pos++;
        } else if(edits[i].type == EDIT_TYPE_REF_GAP) {
            // insert on read
            read.push_back(edits[i].qchr);
            read_pos++;
        } else if(edits[i].type == EDIT_TYPE_MM) {
            assert_eq(ref[ref_pos], edits[i].chr);
            read.push_back(edits[i].qchr);

            read_pos++;
            ref_pos++;
        } else {
            assert(false);
        }
    }

    for(; read_pos < read_len; read_pos++, ref_pos++) {
        read.push_back(ref[ref_pos]);
    }

    return read;
}

template<typename index_t>
static bool compareRepeatCoordByJoinedOff(const RepeatCoord<index_t>& a, const RepeatCoord<index_t>& b)
{
	return a.joinedOff < b.joinedOff;
}


template<typename TStr>
string getString(const TStr& ref, TIndexOffU start, size_t len)
{
    ASSERT_ONLY(const size_t ref_len = ref.length());
    assert_leq(start + len, ref_len);

    string s;
	for(size_t i = 0; i < len; i++) {
		char nt = "ACGT"[ref[start + i]];
		s.push_back(nt);
	}

	return s;
}

template<typename TStr>
void getString(const TStr& ref, TIndexOffU start, size_t len, string& s)
{
    s.clear();
    ASSERT_ONLY(const size_t ref_len = ref.length());
    assert_leq(start + len, ref_len);
    for(size_t i = 0; i < len; i++) {
        char nt = "ACGT"[ref[start + i]];
        s.push_back(nt);
    }
}


template<typename TStr>
inline uint8_t getSequenceBase(const TStr& ref, TIndexOffU pos)
{
    assert_lt(pos, ref.length());
    return (uint8_t)ref[pos];
}


template<typename TStr>
void dump_tstr(const TStr& s)
{
	static int print_width = 60;

	size_t s_len = s.length();

	for(size_t i = 0; i < s_len; i += print_width) {
		string buf;
		for(size_t j = 0; (j < print_width) && (i + j < s_len); j++) {
			buf.push_back("ACGTN"[s[i + j]]);
		}
		cerr << buf << endl;
	}
	cerr << endl;
}


size_t getMaxMatchLen(const EList<Edit>& edits, const size_t read_len)
{
    size_t max_length = 0;
    uint32_t last_edit_pos = 0;
    uint32_t len = 0;

    if (edits.size() == 0) {
        // no edits --> exact match
        return read_len;
    }

    for(size_t i = 0; i < edits.size(); i++) {
        if (last_edit_pos > edits[i].pos) {
            continue;
        }

        len = edits[i].pos - last_edit_pos;
        if(len > max_length) {
            max_length = len;
        }

        last_edit_pos = edits[i].pos + 1;
    }

    if (last_edit_pos < read_len) {
        len = read_len - last_edit_pos;
        if(len > max_length) {
            max_length = len;
        }
    }

    return max_length;
}

int max_index(size_t array[4])
{
    int max_idx = 0;
    
    for(size_t i = 1; i < 4; i++) {
        if(array[max_idx] < array[i]) {
            max_idx = i;
        }
    }
    
    return max_idx;
}

CoordHelper::CoordHelper(TIndexOffU length,
                         TIndexOffU forward_length,
                         const EList<RefRecord>& szs,
                         const EList<string>& ref_names) :
length_(length),
forward_length_(forward_length),
szs_(szs),
ref_namelines_(ref_names)
{
    // build ref_names_ from ref_namelines_
    buildNames();
    buildJoinedFragment();
}

CoordHelper::~CoordHelper()
{
}

void CoordHelper::buildNames()
{
    ref_names_.resize(ref_namelines_.size());
    for(size_t i = 0; i < ref_namelines_.size(); i++) {
        string& nameline = ref_namelines_[i];
        
        for(size_t j = 0; j < nameline.length(); j++) {
            char n = nameline[j];
            if(n == ' ') {
                break;
            }
            ref_names_[i].push_back(n);
        }
    }
}

int CoordHelper::mapJoinedOffToSeq(TIndexOffU joinedOff)
{
    /* search from cached_list */
    if(num_cached_ > 0) {
        for(size_t i = 0; i < num_cached_; i++) {
            Fragments *frag = &cached_[i];
            if(frag->contain(joinedOff)) {
                return frag->frag_id;
            }
        }
        /* fall through */
    }
    
    /* search list */
    int top = 0;
    int bot = fraglist_.size() - 1;
    int pos = 0;
    
    Fragments *frag = &fraglist_[pos];
    while((bot - top) > 1) {
        pos = top + ((bot - top) >> 1);
        frag = &fraglist_[pos];
        
        if (joinedOff < frag->joinedOff) {
            bot = pos;
        } else {
            top = pos;
        }
    }
    
    frag = &fraglist_[top];
    if (frag->contain(joinedOff)) {
        // update cache
        if (num_cached_ < CACHE_SIZE_JOINEDFRG) {
            cached_[num_cached_] = *frag;
            num_cached_++;
        } else {
            cached_[victim_] = *frag;
            victim_++; // round-robin
            victim_ %= CACHE_SIZE_JOINEDFRG;
        }
        
        return top;
    }
    
    return -1;
}

int CoordHelper::getGenomeCoord(TIndexOffU joinedOff,
                                  string& chr_name,
                                  TIndexOffU& pos_in_chr)
{
    int seq_id = mapJoinedOffToSeq(joinedOff);
    if (seq_id < 0) {
        return -1;
    }
    
    Fragments *frag = &fraglist_[seq_id];
    TIndexOffU offset = joinedOff - frag->joinedOff;
    
    pos_in_chr = frag->seqOff + offset;
    chr_name = ref_names_[frag->seq_id];
    
    return 0;
}

void CoordHelper::buildJoinedFragment()
{
    TIndexOffU acc_joinedOff = 0;
    TIndexOffU acc_seqOff = 0;
    TIndexOffU seq_id = 0;
    TIndexOffU frag_id = 0;
    for(size_t i = 0; i < szs_.size(); i++) {
        const RefRecord& rec = szs_[i];
        if(rec.len == 0) {
            continue;
        }
        if (rec.first) {
            acc_seqOff = 0;
            seq_id++;
        }
        
        fraglist_.expand();
        
        fraglist_.back().joinedOff = acc_joinedOff;
        fraglist_.back().length = rec.len;
        
        fraglist_.back().frag_id = frag_id++;
        fraglist_.back().seq_id = seq_id - 1;
        fraglist_.back().first = rec.first;
        
        acc_seqOff += rec.off;
        fraglist_.back().seqOff = acc_seqOff;
        
        acc_joinedOff += rec.len;
        acc_seqOff += rec.len;
    }
    
    // Add Last Fragment(empty)
    fraglist_.expand();
    fraglist_.back().joinedOff = acc_joinedOff;
    fraglist_.back().length = 0;
    fraglist_.back().seqOff = acc_seqOff;
    fraglist_.back().first = false;
    fraglist_.back().frag_id = frag_id;
    fraglist_.back().seq_id = seq_id;
}

TIndexOffU CoordHelper::getEnd(TIndexOffU e) {
    assert_lt(e, length_)
    
    TIndexOffU end = 0;
    if(e < forward_length_) {
        int frag_id = mapJoinedOffToSeq(e);
        assert_geq(frag_id, 0);
        end = fraglist_[frag_id].joinedOff + fraglist_[frag_id].length;
    } else {
        // ReverseComplement
        // a, b are positions w.r.t reverse complement string.
        // fragment map is based on forward string
        int frag_id = mapJoinedOffToSeq(length_ - e - 1);
        assert_geq(frag_id, 0);
        end = length_ - fraglist_[frag_id].joinedOff;
    }
    
    assert_leq(end, length_);
    return end;
}

TIndexOffU CoordHelper::getStart(TIndexOffU e) {
    assert_lt(e, length_)
    
    TIndexOffU start = 0;
    if(e < forward_length_) {
        int frag_id = mapJoinedOffToSeq(e);
        assert_geq(frag_id, 0);
        start = fraglist_[frag_id].joinedOff;
    } else {
        // ReverseComplement
        // a, b are positions w.r.t reverse complement string.
        // fragment map is based on forward string
        int frag_id = mapJoinedOffToSeq(length_ - e - 1);
        assert_geq(frag_id, 0);
        start = length_ - (fraglist_[frag_id].joinedOff + fraglist_[frag_id].length);
    }
    
    assert_leq(start, length_);
    return start;
}

template<typename TStr>
void SeedExt::getExtendedSeedSequence(const TStr& s,
                                      string& seq) const
{
    seq.clear();
    TIndexOffU prev_end = orig_pos.first;
    for(size_t j = 0; j < left_gaps.size(); j++) {
        TIndexOffU curr_end = orig_pos.first - left_gaps[j].first;
        assert_leq(curr_end, prev_end);
        if(curr_end < prev_end) {
            seq = getString(s, curr_end, prev_end - curr_end) + seq;
        }
        int gap_len = left_gaps[j].second;
        assert_neq(gap_len, 0);
        if(gap_len > 0) { // deletion
            string gap_str(gap_len, '-');
            seq = gap_str + seq;
        } else {
            curr_end += gap_len;
        }
        prev_end = curr_end;
    }
    assert_leq(pos.first, prev_end);
    if(pos.first < prev_end) {
        seq = getString(s, pos.first, prev_end - pos.first) + seq;
    }
    if(orig_pos.second - orig_pos.first > 0)
        seq += getString(s, orig_pos.first, orig_pos.second - orig_pos.first);
    TIndexOffU prev_begin = orig_pos.second;
    for(size_t j = 0; j < right_gaps.size(); j++) {
        TIndexOffU curr_begin = orig_pos.second + right_gaps[j].first;
        assert_leq(prev_begin, curr_begin);
        if(prev_begin < curr_begin) {
            seq += getString(s, prev_begin, curr_begin - prev_begin);
        }
        int gap_len = right_gaps[j].second;
        assert_neq(gap_len, 0);
        if(gap_len > 0) { // deletion
            string gap_str(gap_len, '-');
            seq += gap_str;
        } else {
            curr_begin -= gap_len;
        }
        prev_begin = curr_begin;
    }
    assert_leq(prev_begin, pos.second);
    if(prev_begin < pos.second) {
        seq += getString(s, prev_begin, pos.second - prev_begin);
    }
}

SeedSNP *lookup_add_SNP(EList<SeedSNP *>& repeat_snps, SeedSNP& snp)
{
    for(size_t i = 0; i < repeat_snps.size(); i++) {
        if(*repeat_snps[i] == snp) {
            return repeat_snps[i];
        }
    }

    repeat_snps.expand();
    repeat_snps.back() = new SeedSNP();
    *(repeat_snps.back()) = snp;
    return repeat_snps.back();
}

template<typename TStr>
void SeedExt::generateSNPs(const TStr &s, const string& consensus, EList<SeedSNP *>& repeat_snps)
{
    EList<pair<TIndexOffU, int> > gaps;

    // Merge Gaps
    {
        TIndexOffU left_ext_len = getLeftExtLength();
        for(size_t g = 0; g < left_gaps.size(); g++) {
            gaps.expand();
            gaps.back().first = left_ext_len - left_gaps[g].first + left_gaps[g].second;
            gaps.back().second = left_gaps[g].second;
        }
        TIndexOffU right_base = orig_pos.second - pos.first;
        for(size_t g = 0; g < right_gaps.size(); g++) {
            gaps.expand();
            gaps.back().first = right_base + right_gaps[g].first;
            gaps.back().second = right_gaps[g].second;
        }
        gaps.sort();
    }

    //
    string con_ref;
    string seq_read;
    TIndexOffU prev_con_pos = consensus_pos.first;
    TIndexOffU prev_seq_pos = pos.first;

    for(size_t g = 0; g < gaps.size(); g++) {
        TIndexOffU curr_seq_pos = pos.first + gaps[g].first;
        TIndexOffU curr_con_pos = prev_con_pos + (curr_seq_pos - prev_seq_pos);

        con_ref = consensus.substr(prev_con_pos, curr_con_pos - prev_con_pos);
        seq_read = getString(s, prev_seq_pos, curr_seq_pos - prev_seq_pos);
        for(size_t l = 0; l < con_ref.length(); l++) {
            if(seq_read[l] != con_ref[l]) {
                // Single
                SeedSNP snp;
                snp.init(EDIT_TYPE_MM, prev_con_pos + l, 1, seq_read[l]);

                // lookup repeat_snp
                snps.expand();
                snps.back() = lookup_add_SNP(repeat_snps, snp);
            }
        }

        int gap_len = gaps[g].second;
        assert_neq(gap_len, 0);
        if(gap_len > 0) {
            // Deletion (gap on read)

            string gap_str(gap_len, '-');

            SeedSNP snp;
            snp.init(EDIT_TYPE_READ_GAP, curr_con_pos, gap_len, consensus.substr(curr_con_pos, gap_len));

            snps.expand();
            snps.back() = lookup_add_SNP(repeat_snps, snp);

            curr_con_pos += gap_len;

        } else {
            // Insertion (gap on reference)
            string gap_str(abs(gap_len), '-');

            SeedSNP snp;
            snp.init(EDIT_TYPE_REF_GAP, curr_con_pos, abs(gap_len), getString(s, curr_seq_pos, abs(gap_len)));

            snps.expand();
            snps.back() = lookup_add_SNP(repeat_snps, snp);

            curr_seq_pos += abs(gap_len);
        }

        prev_con_pos = curr_con_pos;
        prev_seq_pos = curr_seq_pos;
    }

    assert_eq(consensus_pos.second - prev_con_pos, pos.second - prev_seq_pos);
    con_ref = consensus.substr(prev_con_pos, consensus_pos.second - prev_con_pos);
    seq_read = getString(s, prev_seq_pos, pos.second - prev_seq_pos);

    for(size_t l = 0; l < con_ref.length(); l++) {
        if(seq_read[l] != con_ref[l]) {
            // Single
            
            SeedSNP snp;
            snp.init(EDIT_TYPE_MM, prev_con_pos + l, 1, seq_read[l]);

            snps.expand();
            snps.back() = lookup_add_SNP(repeat_snps, snp);
        }
    }
}

bool seedCmp(const SeedExt& s, const SeedExt& s2)
{
    if(s.getLength() != s2.getLength())
        return s.getLength() > s2.getLength();
    if(s.pos.first != s2.pos.first)
        return s.pos.first < s2.pos.first;
    if(s.pos.second != s2.pos.second)
        return s.pos.second < s2.pos.second;
    return false;
}

template<typename TStr>
void RB_Repeat::init(const RepeatParameter& rp,
                     const TStr& s,
                     CoordHelper& coordHelper,
                     const RB_SubSA& subSA,
                     const RB_RepeatBase& repeatBase)
{
    consensus_ = repeatBase.seq;
    assert_geq(consensus_.length(), rp.min_repeat_len);
    assert_eq(consensus_.length(), rp.min_repeat_len + repeatBase.nodes.size() - 1);
    
    EList<TIndexOffU> positions, next_positions;
    
    const EList<TIndexOffU>& repeat_index = subSA.getRepeatIndex();
    seeds_.clear();
    for(size_t n = 0; n < repeatBase.nodes.size(); n++) {
        TIndexOffU idx = repeatBase.nodes[n];
        TIndexOffU saBegin = repeat_index[idx];
        TIndexOffU saEnd = (idx + 1 < repeat_index.size() ? repeat_index[idx+1] : subSA.size());
        
        next_positions.clear();
        for(size_t sa = saBegin; sa < saEnd; sa++) {
            TIndexOffU left = subSA.get(sa);
            TIndexOffU right = left + subSA.seed_len();
            next_positions.push_back(left);
            
            if(left > 0) {
                size_t idx = positions.bsearchLoBound(left - 1);
                if(idx < positions.size() && positions[idx] == left - 1)
                    continue;
            }
            
            seeds_.expand();
            SeedExt& seed = seeds_.back();
            seed.reset();
            seed.orig_pos = pair<TIndexOffU, TIndexOffU>(left, right);
            seed.pos = seed.orig_pos;
            seed.consensus_pos.first = n;
            seed.consensus_pos.second = n + rp.min_repeat_len;
            seed.bound = pair<TIndexOffU, TIndexOffU>(coordHelper.getStart(left), coordHelper.getEnd(left));
            
#ifndef NDEBUG
            string tmp_str = getString(s, seed.pos.first, seed.pos.second - seed.pos.first);
            assert_eq(tmp_str, consensus_.substr(n, seed.pos.second - seed.pos.first));
#endif

            for(size_t p = seed.consensus_pos.second; p < consensus_.length(); p++) {
                size_t pos = seed.pos.second;
                if(pos >= seed.bound.second)
                    break;
                int nt = getSequenceBase(s, pos);
                if("ACGT"[nt] != consensus_[p])
                    break;
                seed.pos.second++;
                seed.consensus_pos.second++;
            }
        }
        positions = next_positions;
    }
    
    internal_update();
}

template<typename TStr>
void RB_Repeat::extendConsensus(const RepeatParameter& rp,
                                const TStr& s)
{
    size_t remain = seeds_.size();
    EList<string> left_consensuses, right_consensuses, empty_consensuses;
    EList<size_t> ed_seed_nums;

    const TIndexOffU default_max_ext_len = (rp.max_repeat_len - rp.seed_len) / 2;
    const TIndexOffU seed_mm = 1;
    string left_ext_consensus = "", right_ext_consensus = "";
    
    empty_consensuses.resize(seed_mm + 1);
    empty_consensuses.fill("");
    
    while(remain >= rp.repeat_count) {
        for(size_t i = 0; i < remain; i++) {
            seeds_[i].done = false;
            seeds_[i].curr_ext_len = 0;
        }
        
        // extend seeds in left or right direction
        TIndexOffU max_ext_len = min(default_max_ext_len, (TIndexOffU)(rp.max_repeat_len - consensus_.length()));
        get_consensus_seq(s,
                          seeds_,
                          0,      // seed begin
                          remain, // seed end
                          max_ext_len,
                          max_ext_len,
                          seed_mm,
                          rp,
                          ed_seed_nums,
                          &left_consensuses,
                          &right_consensuses);
        
        size_t allowed_seed_mm = 0;
        left_ext_consensus.clear();
        right_ext_consensus.clear();
        for(int i = 0 /* (int)seed_mm */; i >= 0; i--) {
            size_t extlen = (ed_seed_nums[i] < rp.repeat_count ? 0 : (left_consensuses[i].length() + right_consensuses[i].length()));
            // if(extlen <= 0 || extlen < max_ext_len * i / seed_mm)
            //if(extlen < max_ext_len * 2)
            //    continue;
            if(extlen <= 0)
                continue;
            
            left_ext_consensus = left_consensuses[i];
            right_ext_consensus = right_consensuses[i];
            allowed_seed_mm = (size_t)i;
            assert_gt(left_ext_consensus.length(), 0);
            assert_gt(right_ext_consensus.length(), 0);
            break;
        }
        size_t num_passed_seeds = 0;
        if(left_ext_consensus.length() > 0 && right_ext_consensus.length()) {
            consensus_ = reverse(left_ext_consensus) + consensus_;
            consensus_ += right_ext_consensus;
            
#if 0
            calc_edit_dist(s,
                           seeds_,
                           0,
                           remain,
                           left ? consensuses : empty_consensuses,
                           left ? empty_consensuses : consensuses,
                           allowed_seed_mm);
#endif
            
            // update seeds
            for(size_t i = 0; i < seeds_.size(); i++) {
                SeedExt& seed = seeds_[i];
                
                if(i < remain) {
                    if(seed.ed <= allowed_seed_mm) {
                        num_passed_seeds++;
                        seed.done = true;
                        seed.total_ed += seed.ed;
                        seed.pos.first -= left_ext_consensus.length();
                        seed.pos.second += right_ext_consensus.length();
                        seed.consensus_pos.first = 0;
                        seed.consensus_pos.second = consensus_.length();
                        
                        if(seed.left_gaps.size() > 0 &&
                           seed.left_gaps.back().first >= seed.orig_pos.first - seed.pos.first) {
                            int gap_len = seed.left_gaps.back().second;
                            seed.pos.first += gap_len;
                        }
                        if(seed.right_gaps.size() > 0 &&
                           seed.right_gaps.back().first >= seed.pos.second - seed.orig_pos.second) {
                            int gap_len = seed.right_gaps.back().second;
                            seed.pos.second -= gap_len;
                        }
                    }
                }
            }
            
            // move up "done" seeds
            size_t j = 0;
            for(size_t i = 0; i < remain; i++) {
                if(!seeds_[i].done) continue;
                assert_geq(i, j);
                if(i > j) {
                    SeedExt temp = seeds_[j];
                    seeds_[j] = seeds_[i];
                    seeds_[i] = temp;
                    // Find next "undone" seed
                    j++;
                    while(j < i && seeds_[j].done) {
                        j++;
                    }
                    assert(j < remain && !seeds_[j].done);
                } else {
                    j = i + 1;
                }
            }
        }

        remain = num_passed_seeds;
        break;
    } // while(remain >= rp.repeat_count)
    
#ifndef NDEBUG
    // make sure seed positions are unique
    EList<pair<TIndexOffU, TIndexOffU> > seed_poses;
    for(size_t i = 0; i < seeds_.size(); i++) {
        seed_poses.expand();
        seed_poses.back() = seeds_[i].orig_pos;
    }
    seed_poses.sort();
    for(size_t i = 0; i + 1 < seed_poses.size(); i++) {
        if(seed_poses[i].first == seed_poses[i+1].first) {
            assert_lt(seed_poses[i].second, seed_poses[i+1].second);
        } else {
            assert_lt(seed_poses[i].first, seed_poses[i+1].first);
        }
    }
    
    for(size_t i = 0; i < seeds_.size(); i++) {
        assert(seeds_[i].valid());
    }
#endif
    
    RB_Repeat::internal_update();
    
    // check repeats within a repeat sequence
    self_repeat_ = false;
    for(size_t i = 0; i + 1 < seed_ranges_.size(); i++) {
        const RB_AlleleCoord& range = seed_ranges_[i];
        for(size_t j = i + 1; j < seed_ranges_.size(); j++) {
            RB_AlleleCoord& range2 = seed_ranges_[j];
            if(range.right <= range2.left)
                break;
            self_repeat_  = true;
        }
    }
}

template<typename TStr>
void RB_Repeat::getNextRepeat(const RepeatParameter& rp,
                              const TStr& s,
                              RB_Repeat& o)
{
    o.reset();
    
    size_t i = 0;
    for(i = 0; i < seeds_.size() && seeds_[i].done; i++);
    if(i >= seeds_.size())
        return;
    
    for(size_t j = i; j < seeds_.size(); j++) {
        assert(!seeds_[j].done);
        o.seeds_.expand();
        o.seeds_.back() = seeds_[j];
    
    }
    seeds_.resizeExact(i);
    internal_update();
    
    assert_gt(o.seeds_.size(), 0);
    o.consensus_ = getString(s, o.seeds_[0].orig_pos.first, o.seeds_[0].orig_pos.second - o.seeds_[0].orig_pos.first);
    for(size_t j = 0; j < o.seeds_.size(); j++) {
        o.seeds_[j].pos = o.seeds_[j].orig_pos;
        o.seeds_[j].left_gaps.clear();
        o.seeds_[j].right_gaps.clear();
        o.seeds_[j].ed = o.seeds_[j].total_ed = 0;
    }
}

void RB_Repeat::internal_update()
{
    sort(seeds_.begin(), seeds_.end(), seedCmp);
    
    seed_ranges_.resizeExact(seeds_.size());
    for(size_t i = 0; i < seeds_.size(); i++) {
        seed_ranges_[i].left = seeds_[i].pos.first;
        seed_ranges_[i].right = seeds_[i].pos.second;
        seed_ranges_[i].idx = i;
    }
    seed_ranges_.sort();
    
    // check repeats within a repeat sequence
    size_t remove_count = 0;
    for(size_t i = 0; i + 1 < seed_ranges_.size(); i++) {
        const RB_AlleleCoord& range = seed_ranges_[i];
        if(range.left == numeric_limits<size_t>::max())
            continue;
        for(size_t j = i + 1; j < seed_ranges_.size(); j++) {
            RB_AlleleCoord& range2 = seed_ranges_[j];
            if(range2.left == numeric_limits<size_t>::max())
                continue;
            if(range.right <= range2.left)
                break;
            if(range.right >= range2.right) {
                range2.left = numeric_limits<size_t>::max();
                remove_count++;
            }
        }
    }
    
    if(remove_count <= 0)
        return;
    
    for(size_t i = 0; i < seed_ranges_.size(); i++) {
        if(seed_ranges_[i].left == numeric_limits<size_t>::max())
            seeds_[seed_ranges_[i].idx].reset();
    }
    sort(seeds_.begin(), seeds_.end(), seedCmp);
    
#ifndef NDEBUG
    for(size_t i = 0; i < seeds_.size(); i++) {
        if(i < seeds_.size() - remove_count) {
            assert_lt(seeds_[i].pos.first, seeds_[i].pos.second);
        } else {
            assert_eq(seeds_[i].pos.first, seeds_[i].pos.second);
        }
    }
#endif
    
    seeds_.resize(seeds_.size() - remove_count);
    seed_ranges_.resize(seeds_.size());
    for(size_t i = 0; i < seeds_.size(); i++) {
        seed_ranges_[i].left = seeds_[i].pos.first;
        seed_ranges_[i].right = seeds_[i].pos.second;
        seed_ranges_[i].idx = i;
    }
    seed_ranges_.sort();
}

bool RB_Repeat::contain(TIndexOffU left, TIndexOffU right) const
{
    size_t l = 0, r = seed_ranges_.size();
    while(l < r) {
        size_t m = (l + r) / 2;
        const RB_AlleleCoord& coord = seed_ranges_[m];
        if(right <= coord.left) {
            r = m;
        } else if(left >= coord.right) {
            l = m + 1;
        } else {
            return coord.left <= left && right <= coord.right;
        }
    }
    
    return false;
}

template<typename TStr>
void RB_Repeat::saveSeedExtension(const RepeatParameter& rp,
                                  const TStr& s,
                                  CoordHelper& coordHelper,
                                  TIndexOffU grp_id,
                                  ostream& fp,
                                  size_t& total_repeat_seq_len,
                                  size_t& total_allele_seq_len) const
{
    // apply color, which is compatible with linux commands such as cat and less -r
#if 1
    const string red = "\033[31m", reset = "\033[0m", resetline = "\033[4m", redline = "\033[4;31m";
#else
    const string red = "", reset = "", resetline = "", redline = "";
#endif
    bool show_exterior_seq = true;
    const size_t max_show_seq_len = 700;
    
    size_t total_count = 0;
    for(size_t i = 0; i < seeds_.size(); i++) {
        const SeedExt& seed = seeds_[i];
        size_t ext_len = seed.getLength();
        if(ext_len < rp.min_repeat_len) continue;
        total_allele_seq_len += ext_len;
        total_count++;
        bool sense_strand = seed.pos.first < coordHelper.forward_length();
        
        fp << setw(6) << repeat_id_ << "  " << setw(5) << seeds_.size();
        fp << "  " << setw(4) << i;
        fp << "  " << setw(4) << ext_len;
        fp << "  " << (sense_strand ? '+' : '-');
        fp << "  " << setw(10) << seed.pos.first << "  " << setw(10) << seed.pos.second;
        fp << "  " << setw(10) << seed.orig_pos.first << "  " << setw(10) << seed.orig_pos.second;
        
        string chr_name;
        TIndexOffU pos_in_chr;
        if(sense_strand) {
            coordHelper.getGenomeCoord(seed.pos.first, chr_name, pos_in_chr);
        } else {
            coordHelper.getGenomeCoord(s.length() - seed.pos.first - (seed.pos.second - seed.pos.first), chr_name, pos_in_chr);
        }
        fp << "  " << setw(5) << chr_name << ":" << setw(10) << std::left << pos_in_chr << std::right;
        
        if(sense_strand) {
            coordHelper.getGenomeCoord(seed.pos.second, chr_name, pos_in_chr);
        } else {
            coordHelper.getGenomeCoord(s.length() - seed.pos.second - (seed.pos.second - seed.pos.first), chr_name, pos_in_chr);
        }
        fp << "  " << setw(5) << chr_name << ":" << setw(10) << std::left << pos_in_chr << std::right;
        
        if(!seed.aligned) {
            fp << endl;
            continue;
        }
        
        string deststr = "";
        seed.getExtendedSeedSequence(s, deststr);
        
        // add exterior sequences
        if(seed.consensus_pos.first > 0) {
            TIndexOffU seq_pos, seq_len;
            if(seed.pos.first >= seed.consensus_pos.first) {
                seq_pos = seed.pos.first - seed.consensus_pos.first;
                seq_len = seed.consensus_pos.first;
            } else {
                seq_pos = 0;
                seq_len = seed.pos.first;
            }
            deststr = getString(s, seq_pos, seq_len) + deststr;
            if(seq_len < seed.consensus_pos.first) {
                deststr = string(seed.consensus_pos.first - seq_len, 'N') + deststr;
            }
        }
        if(seed.consensus_pos.second < consensus_.length()) {
            deststr += getString(s, seed.pos.second, consensus_.length() - seed.consensus_pos.second);
        }
        
        assert_eq(consensus_.length(), deststr.length());
        fp << "  ";
        
        // print sequence w.r.t. the current group
        for(size_t j = 0; j < min(consensus_.length(), max_show_seq_len); j++) {
            bool outside = j < seed.consensus_pos.first || j >= seed.consensus_pos.second;
            bool different = (consensus_[j] != deststr[j]);
            if(outside) {
                if(show_exterior_seq) {
                    if(different) {
                        fp << redline;
                        fp << deststr[j];
                        fp << reset;
                    } else {
                        fp << resetline;
                        fp << deststr[j];
                        fp << reset;
                    }
                } else {
                    if(j < seed.consensus_pos.first)
                        fp << " ";
                }
            } else {
                if(different) fp << red;
                fp << deststr[j];
                if(different) fp << reset;
            }
        }
        
#if 0
        fp << "\t";
        for(size_t ei = 0; ei < seed.edits.size(); ei++) {
            const Edit& edit = seed.edits[ei];
            if(ei > 0) fp << ",";
            fp << edit;
            if (edit.snpID != std::numeric_limits<uint32_t>::max()) {
                fp << "@" << edit.snpID;
            }
        }
#endif
        
        fp << endl;
    }
    
    if(total_count > 0) fp << setw(5) << total_count << endl << endl;
    total_repeat_seq_len += consensus_.length();
}

template<typename TStr>
size_t calc_edit_dist(const TStr& s,
                      const SeedExt& seed,
                      const SeedExt& seed2,
                      size_t left_ext,
                      size_t right_ext,
                      size_t max_ed)
{
    if(seed.bound.first + left_ext > seed.pos.first ||
       seed.pos.second + right_ext > seed.bound.second ||
       seed2.bound.first + left_ext > seed2.pos.first ||
       seed2.pos.second + right_ext > seed2.bound.second)
        return max_ed + 1;
    
    size_t ed = 0;
    for(size_t i = 0; i < left_ext; i++) {
        int ch = getSequenceBase(s, seed.pos.first - i - 1);
        assert_range(0, 3, ch);
        int ch2 = getSequenceBase(s, seed2.pos.first - i - 1);
        assert_range(0, 3, ch2);
        if(ch != ch2) ed++;
        if(ed > max_ed)
            return ed;
        
    }
    for(size_t i = 0; i < right_ext; i++) {
        int ch = getSequenceBase(s, seed.pos.second + i);
        assert_range(0, 3, ch);
        int ch2 = getSequenceBase(s, seed2.pos.second + i);
        assert_range(0, 3, ch2);
        if(ch != ch2) ed++;
        if(ed > max_ed)
            return ed;
    }
    
    return ed;
}

size_t extract_kmer(const string& seq, size_t offset, size_t k = 5)
{
    assert_leq(offset + k, seq.length());
    size_t kmer = 0;
    for(size_t i = 0; i < k; i++) {
        kmer = (kmer << 2) | asc2dna[seq[offset + i]];
    }
    return kmer;
}

size_t next_kmer(size_t kmer, char base, size_t k = 5)
{
    kmer &= ((1 << ((k-1)*2))) - 1;
    kmer = (kmer << 2) | asc2dna[base];
    return kmer;
}

void build_kmer_table(const string& consensus, EList<pair<size_t, size_t> >& kmer_table, size_t k = 5)
{
    kmer_table.clear();
    if(consensus.length() < k)
        return;
    size_t kmer = 0;
    for(size_t i = 0; i + k <= consensus.length(); i++) {
        if(i == 0) {
            kmer = extract_kmer(consensus, i, k);
        } else {
            kmer = next_kmer(kmer, consensus[i+k-1], k);
        }
        kmer_table.expand();
        kmer_table.back().first = kmer;
        kmer_table.back().second = i;
    }
    kmer_table.sort();
}

void find_gap_pos(const string& s,
                  const string& s2,
                  EList<size_t>& ed,
                  EList<size_t>& ed2,
                  bool del,
                  size_t gap_len,
                  size_t max_mm,
                  size_t& gap_pos,
                  size_t& mm,
                  bool debug = false)
{
    assert_eq(s.length(), s2.length());
    size_t seq_len = s.length();
    ed.resizeExact(seq_len); ed.fill(max_mm + 1);
    ed2.resizeExact(seq_len); ed2.fill(max_mm + 1);
    
    // from left to right
    for(size_t i = 0; i < seq_len; i++) {
        size_t add = (s[i] == s2[i] ? 0 : 1);
        if(i == 0) ed[i] = add;
        else       ed[i] = ed[i-1] + add;
        if(ed[i] >= max_mm + 1) break;
    }
    
    // from right to left
    size_t s_sub  = del ? 0 : gap_len;
    size_t s2_sub = del ? gap_len : 0;
    for(int i = seq_len - 1; i >= gap_len; i--) {
        size_t add = (s[i - s_sub] == s2[i - s2_sub] ? 0 : 1);
        if(i == seq_len - 1) ed2[i] = add;
        else                 ed2[i] = ed2[i+1] + add;
        if(ed2[i] > max_mm) break;
    }
    
    if(debug) {
        cout << s << endl << s2 << endl;
        for(size_t i = 0; i < ed.size(); i++) {
            cout << (ed[i] % 10);
        }
        cout << endl;
        for(size_t i = 0; i < ed2.size(); i++) {
            cout << (ed2[i] % 10);
        }
        cout << endl;
    }
    
    size_t min_mm = ed2[gap_len];
    int min_mm_i = -1;
    assert_eq(ed.size(), ed2.size());
    for(size_t i = 0; i + gap_len + 1 < ed.size(); i++) {
        if(ed[i] > max_mm) break;
        size_t cur_mm = ed[i] + ed2[i + gap_len + 1];
        if(cur_mm < min_mm) {
            min_mm = cur_mm;
            min_mm_i = i;
        }
    }
    
    mm = min_mm;
    if(mm > max_mm)
        return;
    gap_pos = (size_t)(min_mm_i + 1);
}

void align_with_one_gap(const string& s,
                        const EList<pair<size_t, size_t> >& s_kmer_table,
                        const string& s2,
                        EList<size_t>& counts,
                        size_t max_gap,
                        size_t max_mm,
                        size_t& mm,
                        size_t& gap_pos,
                        int& gap_len,
                        size_t k = 5)
{
    mm = max_mm + 1;
    assert_eq(s.length(), s2.length());
    counts.resizeExact(max_gap * 2 + 1);
    counts.fillZero();
    size_t max_count = 0, max_count_i = 0;
    for(size_t i = 0; i + k <= s2.length(); i += 1) {
        pair<size_t, size_t> kmer(0, 0);
        kmer.first = extract_kmer(s2, i, k);
        size_t lb = s_kmer_table.bsearchLoBound(kmer);
        while(lb < s_kmer_table.size() && kmer.first == s_kmer_table[lb].first) {
            int gap = (int)s_kmer_table[lb].second - (int)i;
            if(gap != 0 && abs(gap) < max_gap) {
                size_t gap_i = (size_t)(gap + max_gap);
                counts[gap_i]++;
                if(counts[gap_i] > max_count) {
                    max_count = counts[gap_i];
                    max_count_i = gap_i;
                }
            }
            lb++;
        }
    }
    
    if(max_count <= 0)
        return;
    
    assert_lt(max_count_i, counts.size());
    int gap = (int)max_count_i - (int)max_gap;
    assert_neq(gap, 0);
    size_t abs_gap = abs(gap);
    bool del = gap > 0;
    
    EList<size_t> ed, ed2;
    find_gap_pos(s,
                 s2,
                 ed,
                 ed2,
                 del,
                 abs_gap,
                 max_mm,
                 gap_pos,
                 mm);
    if(mm > max_mm)
        return;
    
    gap_len = del ? abs_gap : -abs_gap;
    
#ifndef NDEBUG
    assert_leq(mm, max_mm);
    string ds = s, ds2 = s2;
    size_t debug_mm = 0;
    if(del) ds.erase(gap_pos, abs_gap);
    else    ds2.erase(gap_pos, abs_gap);
    
    for(size_t i = 0; i < min(ds.length(), ds2.length()); i++) {
        if(ds[i] != ds2[i])
            debug_mm++;
    }
    assert_eq(debug_mm, mm);
#endif
}

template<typename TStr>
void calc_edit_dist(const TStr& s,
                    EList<SeedExt>& seeds,
                    size_t sb,
                    size_t se,
                    const EList<string>& left_consensuses,
                    const EList<string>& right_consensuses,
                    uint32_t max_ed)
{
    string left_consensus = left_consensuses[max_ed];
    string right_consensus = right_consensuses[max_ed];
    
    EList<pair<size_t, size_t> > left_kmer_table, right_kmer_table;
    build_kmer_table(left_consensus, left_kmer_table);
    build_kmer_table(right_consensus, right_kmer_table);

    string left_seq, right_seq;
    const size_t max_gap = 10;
    EList<size_t> counts;
    
    size_t left_ext = left_consensus.length();
    size_t right_ext = right_consensus.length();
    for(size_t i = sb; i < se; i++) {
        SeedExt& seed = seeds[i];
        if(seed.bound.first + left_ext > seed.pos.first ||
           seed.pos.second + right_ext > seed.bound.second) {
            seed.ed = max_ed + 1;
            continue;
        }
        
        size_t left_ed = 0;
        if(left_ext > 0) {
            getString(s, seed.pos.first - left_ext, left_ext, left_seq);
            reverse(left_seq.begin(), left_seq.end());
            for(size_t j = 0; j < left_ext; j++) {
                if(left_seq[j] != left_consensus[j]) left_ed++;
                if(left_ed <= max_ed && j < left_consensuses[left_ed].length()) {
                    seed.curr_ext_len = j + 1;
                }
            }
            
            if(left_ed > max_ed) {
                size_t gap_pos = 0;
                int gap_len = 0;
                align_with_one_gap(left_consensus,
                                   left_kmer_table,
                                   left_seq,
                                   counts,
                                   min(left_ed - max_ed, max_gap),
                                   max_ed,
                                   left_ed,
                                   gap_pos,
                                   gap_len);
                if(left_ed <= max_ed) {
                    seed.left_gaps.expand();
                    seed.left_gaps.back().first = seed.getLeftExtLength() + gap_pos;
                    seed.left_gaps.back().second = gap_len;
                }
            }
        } else {
            left_seq.clear();
        }
        
        size_t right_ed = 0;
        if(right_ext > 0) {
            getString(s, seed.pos.second, right_ext, right_seq);
            for(size_t j = 0; j < right_ext; j++) {
                if(right_seq[j] != right_consensus[j]) right_ed++;
                if(right_ed <= max_ed && j < right_consensuses[right_ed].length()) {
                    seed.curr_ext_len = j + 1;
                }
            }
            if(right_ed > max_ed) {
                size_t gap_pos = 0;
                int gap_len = 0;
                align_with_one_gap(right_consensus,
                                   right_kmer_table,
                                   right_seq,
                                   counts,
                                   min(right_ed - max_ed, max_gap),
                                   max_ed,
                                   right_ed,
                                   gap_pos,
                                   gap_len);
                if(right_ed <= max_ed) {
                    seed.right_gaps.expand();
                    seed.right_gaps.back().first = seed.getRightExtLength() + gap_pos;
                    seed.right_gaps.back().second = gap_len;
                }
            }
        } else {
            right_seq.clear();
        }
        
        seed.ed = left_ed + right_ed;
    }
}

template<typename TStr>
void RB_Repeat::get_consensus_seq(const TStr& s,
                                  EList<SeedExt>& seeds,
                                  size_t sb,
                                  size_t se,
                                  size_t min_left_ext,
                                  size_t min_right_ext,
                                  size_t max_ed,
                                  const RepeatParameter& rp,
                                  EList<size_t>& ed_seed_nums,
                                  EList<string>* left_consensuses,
                                  EList<string>* right_consensuses) const
{
    assert_lt(sb, se);
    assert_geq(se - sb, rp.seed_count);
    assert_leq(se, seeds.size());
    if(left_consensuses != NULL) {
        left_consensuses->clear(); left_consensuses->resize(max_ed + 1);
        for(size_t i = 0; i < max_ed + 1; i++) {
            (*left_consensuses)[i].clear();
        }
    }
    if(right_consensuses != NULL) {
        right_consensuses->clear(); right_consensuses->resize(max_ed + 1);
        for(size_t i = 0; i < max_ed + 1; i++) {
            (*right_consensuses)[i].clear();
        }
    }
    
    assert_eq(min_left_ext, min_right_ext);
    
    EList<string> seqs;
    seqs.reserveExact(seeds.size());
    size_t max_i = 0, max_count = 0;
    while(min_left_ext > 0) {
        seqs.clear();
        max_i = 0;
        max_count = 0;
        for(size_t i = sb; i < se; i++) {
            const SeedExt& seed = seeds[i];
            if(seeds[i].bound.first + min_left_ext > seed.pos.first ||
               seed.pos.second + min_right_ext > seed.bound.second)
                continue;
            
            seqs.expand();
            if(min_left_ext > 0) {
                getString(s, seed.pos.first - min_left_ext, min_left_ext, seqs.back());
            }
            if(min_right_ext > 0) {
                seqs.back() += getString(s, seed.pos.second, min_right_ext);
            }
        }
        seqs.sort();
        for(size_t i = 0; i + max_count < seqs.size();) {
            size_t count = 1;
            for(size_t j = i + 1; j < seqs.size(); j++) {
                if(seqs[i] == seqs[j]) count++;
                else                   break;
            }
            if(count >= max_count) {
                max_count = count;
                max_i = i;
            }
            i = i + count;
        }
        if(max_count >= rp.seed_count)
            break;
        
        min_left_ext--;
        min_right_ext--;
    }
    
    // update ed
    // extend consensus string
    ed_seed_nums.resize(max_ed + 1); ed_seed_nums.fillZero();
    EList<size_t> next_ed_seed_nums; next_ed_seed_nums.resize(max_ed + 1);
    
    if(max_count < rp.seed_count)
        return;
    
    for(size_t i = sb; i < se; i++) seeds[i].ed = 0;
    size_t seed_ext_len = 0;
    while(seed_ext_len < max(min_left_ext, min_right_ext)) {
        // select base to be used for extension
        uint8_t left_ext_base = 0;
        if(seed_ext_len < min_left_ext) left_ext_base = seqs[max_i][min_left_ext - seed_ext_len - 1];
        uint8_t right_ext_base = 0;
        if(seed_ext_len < min_right_ext) right_ext_base = seqs[max_i][min_left_ext + seed_ext_len];
        
        // estimate extended ed
        next_ed_seed_nums.fillZero();
        for(size_t i = sb; i < se; i++) {
            uint32_t est_ed = seeds[i].ed;
            if(seed_ext_len < min_left_ext) {
                if(seeds[i].bound.first + seed_ext_len + 1 <= seeds[i].pos.first) {
                    TIndexOffU left_pos = seeds[i].pos.first - seed_ext_len - 1;
                    int ch = getSequenceBase(s, left_pos);
                    assert_range(0, 3, ch);
                    if ("ACGT"[ch] != left_ext_base) {
                        est_ed++;
                    }
                } else {
                    est_ed = max_ed + 1;
                }
            }
            
            if(seed_ext_len < min_right_ext) {
                TIndexOffU right_pos = seeds[i].pos.second + seed_ext_len;
                if(right_pos >= seeds[i].bound.second) {
                    est_ed = max_ed + 1;
                } else {
                    int ch = getSequenceBase(s, right_pos);
                    assert_range(0, 3, ch);
                    if ("ACGT"[ch] != right_ext_base) {
                        est_ed++;
                    }
                }
            }
            
            if(est_ed <= max_ed) {
                next_ed_seed_nums[est_ed]++;
            }
        }
        
        for(size_t i = 1; i < next_ed_seed_nums.size(); i++) next_ed_seed_nums[i] += next_ed_seed_nums[i-1];
        if(next_ed_seed_nums[max_ed] < rp.repeat_count) {
            break;
        }
        
        for(size_t i = sb; i < se; i++) {
            if(seed_ext_len < min_left_ext) {
                if(seeds[i].bound.first + seed_ext_len + 1 <= seeds[i].pos.first) {
                    TIndexOffU left_pos = seeds[i].pos.first - seed_ext_len - 1;
                    int ch = getSequenceBase(s, left_pos);
                    assert_range(0, 3, ch);
                    if ("ACGT"[ch] != left_ext_base) {
                        seeds[i].ed++;
                    }
                } else {
                    seeds[i].ed = max_ed + 1;
                }
            }
            
            if(seed_ext_len < min_right_ext) {
                TIndexOffU right_pos = seeds[i].pos.second + seed_ext_len;
                if(right_pos >= seeds[i].bound.second) {
                    seeds[i].ed = max_ed + 1;
                } else {
                    int ch = getSequenceBase(s, right_pos);
                    assert_range(0, 3, ch);
                    if ("ACGT"[ch] != right_ext_base) {
                        seeds[i].ed++;
                    }
                }
            }
        }
        
        for(size_t i = 0; i < next_ed_seed_nums.size(); i++) {
            if(next_ed_seed_nums[i] < rp.repeat_count)
                continue;
            ed_seed_nums[i] = next_ed_seed_nums[i];
            if(seed_ext_len < min_left_ext) {
                assert(left_consensuses != NULL);
                (*left_consensuses)[i] += left_ext_base;
            }
            if(seed_ext_len < min_right_ext) {
                assert(right_consensuses != NULL);
                (*right_consensuses)[i] += right_ext_base;
            }
        }
        
        seed_ext_len++;
    }
}

template<typename TStr>
void RB_RepeatExt::get_consensus_seq(const TStr& s,
                                     EList<SeedExt>& seeds,
                                     size_t sb,
                                     size_t se,
                                     size_t min_left_ext,
                                     size_t min_right_ext,
                                     size_t max_ed,
                                     const RepeatParameter& rp,
                                     EList<size_t>& ed_seed_nums,
                                     EList<string>* left_consensuses,
                                     EList<string>* right_consensuses) const
{
    assert_lt(sb, se);
    assert_leq(se, seeds.size());
    if(left_consensuses != NULL) {
        left_consensuses->clear(); left_consensuses->resize(max_ed + 1);
        for(size_t i = 0; i < max_ed + 1; i++) {
            (*left_consensuses)[i].clear();
        }
    }
    if(right_consensuses != NULL) {
        right_consensuses->clear(); right_consensuses->resize(max_ed + 1);
        for(size_t i = 0; i < max_ed + 1; i++) {
            (*right_consensuses)[i].clear();
        }
    }
    
    // cluster sequences
    EList<size_t> belongto;
    belongto.resizeExact(se - sb);
    for(size_t i = 0; i < belongto.size(); i++) belongto[i] = i;
    
    for(size_t i = 0; i + 1 < belongto.size(); i++) {
        for(size_t j = i + 1; j < belongto.size(); j++) {
            if(belongto[j] != j) continue;
            size_t ed = calc_edit_dist(s,
                                       seeds[sb + i],
                                       seeds[sb + j],
                                       min_left_ext,
                                       min_right_ext,
                                       max_ed + 1);
            if(ed <= max_ed + 1) {
                belongto[j] = belongto[i];
            }
            
        }
    }
    
    // find the maximum group that has the most sequences
    EList<size_t> vote;
    vote.resizeExact(seeds.size());
    vote.fillZero();
    size_t max_group = 0;
    for(size_t i = 0; i < belongto.size(); i++) {
        size_t cur_group = belongto[i];
        assert_lt(cur_group, vote.size());
        vote[cur_group]++;
        if(cur_group != max_group && vote[cur_group] > vote[max_group]) {
            max_group = cur_group;
        }
    }
    
    // reuse vote to store seeds
    EList<size_t>& consensus_group = vote;
    consensus_group.clear();
    for(size_t i = 0; i < belongto.size(); i++) {
        if(belongto[i] == max_group)
            consensus_group.push_back(i);
    }
    
    // update ed
    // extend consensus string
    ed_seed_nums.resize(max_ed + 1); ed_seed_nums.fillZero();
    EList<size_t> next_ed_seed_nums; next_ed_seed_nums.resize(max_ed + 1);
    for(size_t i = sb; i < se; i++) seeds[i].ed = 0;
    size_t seed_ext_len = 0;
    while(seed_ext_len < max(min_left_ext, min_right_ext)) {
        // count base
        size_t l_count[4] = {0, };
        size_t r_count[4] = {0, };
        
        for(size_t i = 0; i < consensus_group.size(); i++) {
            size_t si = sb + consensus_group[i];
            int ch;
            if(seed_ext_len < min_left_ext) {
                if(seeds[i].bound.first + seed_ext_len + 1 <= seeds[i].pos.first) {
                    ch = getSequenceBase(s, seeds[si].pos.first - seed_ext_len - 1);
                    assert_range(0, 3, ch);
                    l_count[ch]++;
                }
            }
            if(seed_ext_len < min_right_ext) {
                if(seeds[i].bound.second + seed_ext_len) {
                    ch = getSequenceBase(s, seeds[si].pos.second + seed_ext_len);
                    assert_range(0, 3, ch);
                    r_count[ch]++;
                }
            }
        }
        
        // select base to be used for extension
        uint8_t left_ext_base = 0;
        if(seed_ext_len < min_left_ext) left_ext_base = max_index(l_count);
        uint8_t right_ext_base = 0;
        if(seed_ext_len < min_right_ext) right_ext_base = max_index(r_count);
        
        // estimate extended ed
        next_ed_seed_nums.fillZero();
        for(size_t i = sb; i < se; i++) {
            uint32_t est_ed = seeds[i].ed;
            if(seed_ext_len < min_left_ext) {
                if(seeds[i].bound.first + seed_ext_len + 1 <= seeds[i].pos.first) {
                    TIndexOffU left_pos = seeds[i].pos.first - seed_ext_len - 1;
                    int ch = getSequenceBase(s, left_pos);
                    assert_range(0, 3, ch);
                    if (ch != left_ext_base) {
                        est_ed++;
                    }
                } else {
                    est_ed = max_ed + 1;
                }
            }
            
            if(seed_ext_len < min_right_ext) {
                TIndexOffU right_pos = seeds[i].pos.second + seed_ext_len;
                if(right_pos >= seeds[i].bound.second) {
                    est_ed = max_ed + 1;
                } else {
                    int ch = getSequenceBase(s, right_pos);
                    assert_range(0, 3, ch);
                    if (ch != right_ext_base) {
                        est_ed++;
                    }
                }
            }
            
            if(est_ed <= max_ed) {
                next_ed_seed_nums[est_ed]++;
            }
        }
        
        for(size_t i = 1; i < next_ed_seed_nums.size(); i++) next_ed_seed_nums[i] += next_ed_seed_nums[i-1];
        if(next_ed_seed_nums[max_ed] < rp.repeat_count) {
            break;
        }
        
        for(size_t i = sb; i < se; i++) {
            if(seed_ext_len < min_left_ext) {
                if(seeds[i].bound.first + seed_ext_len + 1 <= seeds[i].pos.first) {
                    TIndexOffU left_pos = seeds[i].pos.first - seed_ext_len - 1;
                    int ch = getSequenceBase(s, left_pos);
                    assert_range(0, 3, ch);
                    if (ch != left_ext_base) {
                        seeds[i].ed++;
                    }
                } else {
                    seeds[i].ed = max_ed + 1;
                }
            }
            
            if(seed_ext_len < min_right_ext) {
                TIndexOffU right_pos = seeds[i].pos.second + seed_ext_len;
                if(right_pos >= seeds[i].bound.second) {
                    seeds[i].ed = max_ed + 1;
                } else {
                    int ch = getSequenceBase(s, right_pos);
                    assert_range(0, 3, ch);
                    if (ch != right_ext_base) {
                        seeds[i].ed++;
                    }
                }
            }
        }
        
        for(size_t i = 0; i < next_ed_seed_nums.size(); i++) {
            if(next_ed_seed_nums[i] < rp.repeat_count)
                continue;
            ed_seed_nums[i] = next_ed_seed_nums[i];
            if(seed_ext_len < min_left_ext) {
                assert(left_consensuses != NULL);
                (*left_consensuses)[i] += (char)("ACGT"[left_ext_base]);
            }
            if(seed_ext_len < min_right_ext) {
                assert(right_consensuses != NULL);
                (*right_consensuses)[i] += (char)("ACGT"[right_ext_base]);
            }
        }
        
        seed_ext_len++;
    }
}

EList<size_t> RB_Repeat::ca_ed_;
EList<size_t> RB_Repeat::ca_ed2_;
string RB_Repeat::ca_s_;
string RB_Repeat::ca_s2_;

size_t RB_Repeat::seed_merge_tried = 0;
size_t RB_Repeat::seed_merged = 0;

template<typename TStr>
void RB_RepeatExt::extendConsensus(const RepeatParameter& rp,
                                   const TStr& s)
{
    size_t remain = seeds_.size();
    EList<string> consensuses, empty_consensuses;
    EList<size_t> ed_seed_nums;
    bool left = true;
    
    const TIndexOffU default_max_ext_len = 100;
    const TIndexOffU seed_mm = 4;
    string ext_consensus = "";
    
    empty_consensuses.resize(seed_mm + 1);
    empty_consensuses.fill("");
    
    while(remain >= rp.repeat_count) {
        for(size_t i = 0; i < remain; i++) {
            seeds_[i].done = false;
            seeds_[i].curr_ext_len = 0;
        }
        
        // extend seeds in left or right direction
        TIndexOffU max_ext_len = min(default_max_ext_len, (TIndexOffU)(rp.max_repeat_len - consensus_.length()));
        get_consensus_seq(s,
                          seeds_,
                          0,      // seed begin
                          remain, // seed end
                          left ? max_ext_len : 0,
                          left ? 0 : max_ext_len,
                          seed_mm,
                          rp,
                          ed_seed_nums,
                          left ? &consensuses : NULL,
                          left ? NULL : &consensuses);
        
        size_t allowed_seed_mm = 0;
        ext_consensus.clear();
        for(int i = (int)seed_mm; i >= 0; i--) {
            size_t extlen = (ed_seed_nums[i] < rp.repeat_count ? 0 : consensuses[i].length());
            if(extlen <= 0 || extlen < max_ext_len * i / seed_mm)
                continue;
            
            if(i > 0 && consensuses[i].length() <= consensuses[i-1].length() + 5)
                continue;
            
            ext_consensus = consensuses[i];
            allowed_seed_mm = (size_t)i;
            assert(ext_consensus.length() > 0);
            break;
        }
        size_t num_passed_seeds = 0;
        if(ext_consensus.length() > 0) {
            if(left) consensus_ = reverse(ext_consensus) + consensus_;
            else consensus_ += ext_consensus;
            
            calc_edit_dist(s,
                           seeds_,
                           0,
                           remain,
                           left ? consensuses : empty_consensuses,
                           left ? empty_consensuses : consensuses,
                           allowed_seed_mm);
            
            // update seeds
            for(size_t i = 0; i < seeds_.size(); i++) {
                SeedExt& seed = seeds_[i];
                
                if(i < remain) {
                    if(seed.ed <= allowed_seed_mm) {
                        num_passed_seeds++;
                        seed.done = true;
                        seed.total_ed += seed.ed;
                        if(left) {
                            if(seed.left_gaps.size() > 0 &&
                               seed.left_gaps.back().first >= seed.orig_pos.first - seed.pos.first) {
                                int gap_len = seed.left_gaps.back().second;
                                seed.pos.first += gap_len;
                            }
                            seed.pos.first -= ext_consensus.length();
                            seed.consensus_pos.first = 0;
                            seed.consensus_pos.second = consensus_.length();
                        } else {
                            if(seed.right_gaps.size() > 0 &&
                               seed.right_gaps.back().first >= seed.pos.second - seed.orig_pos.second) {
                                int gap_len = seed.right_gaps.back().second;
                                seed.pos.second -= gap_len;
                            }
                            seed.pos.second += ext_consensus.length();
                            seed.consensus_pos.second = consensus_.length();
                        }
                    } else {
                        if(left) {
                            assert_leq(seed.curr_ext_len, ext_consensus.length());
                            TIndexOffU adjust = ext_consensus.length() - seed.curr_ext_len;
                            seed.consensus_pos.first += adjust;
                            seed.consensus_pos.second += ext_consensus.length();
                            assert_leq(seed.curr_ext_len, seed.pos.first);
                            seed.pos.first -= seed.curr_ext_len;
                        } else {
                            assert_leq(seed.curr_ext_len, ext_consensus.length());
                            seed.consensus_pos.second += seed.curr_ext_len;
                            seed.pos.second += seed.curr_ext_len;
                        }
                    }
                } else {
                    if(left) {
                        seed.consensus_pos.first += ext_consensus.length();
                        seed.consensus_pos.second += ext_consensus.length();
                    }
                }
            }
            
            // move up "done" seeds
            size_t j = 0;
            for(size_t i = 0; i < remain; i++) {
                if(!seeds_[i].done) continue;
                assert_geq(i, j);
                if(i > j) {
                    SeedExt temp = seeds_[j];
                    seeds_[j] = seeds_[i];
                    seeds_[i] = temp;
                    // Find next "undone" seed
                    j++;
                    while(j < i && seeds_[j].done) {
                        j++;
                    }
                    assert(j < remain && !seeds_[j].done);
                } else {
                    j = i + 1;
                }
            }
        }
        
        remain = num_passed_seeds;
        if(remain < rp.repeat_count) {
            if(left) {
                left = false;
                remain = seeds_.size();
            }
        }
    } // while(remain >= rp.repeat_count)
    
#ifndef NDEBUG
    // make sure seed positions are unique
    EList<pair<TIndexOffU, TIndexOffU> > seed_poses;
    for(size_t i = 0; i < seeds_.size(); i++) {
        seed_poses.expand();
        seed_poses.back() = seeds_[i].orig_pos;
    }
    seed_poses.sort();
    for(size_t i = 0; i + 1 < seed_poses.size(); i++) {
        if(seed_poses[i].first == seed_poses[i+1].first) {
            assert_lt(seed_poses[i].second, seed_poses[i+1].second);
        } else {
            assert_lt(seed_poses[i].first, seed_poses[i+1].first);
        }
    }
    
    for(size_t i = 0; i < seeds_.size(); i++) {
        assert(seeds_[i].valid());
    }
#endif
    
    internal_update();
    
    // check repeats within a repeat sequence
    self_repeat_ = false;
    for(size_t i = 0; i + 1 < seed_ranges_.size(); i++) {
        const RB_AlleleCoord& range = seed_ranges_[i];
        for(size_t j = i + 1; j < seed_ranges_.size(); j++) {
            RB_AlleleCoord& range2 = seed_ranges_[j];
            if(range.right <= range2.left)
                break;
            self_repeat_  = true;
        }
    }
}

bool RB_Repeat::overlap(const RB_Repeat& o,
                        bool& contain,
                        bool& left,
                        size_t& seed_i,
                        size_t& seed_j,
                        bool debug) const
{
    contain = left = false;
    seed_i = seed_j = 0;
    size_t p = 0, p2 = 0;
    while(p < seed_ranges_.size() && p2 < o.seed_ranges_.size()) {
        const RB_AlleleCoord& range = seed_ranges_[p];
        const SeedExt& seed = seeds_[range.idx];
        Range range_ = seed.getExtendedRange(consensus_.length());
        RB_AlleleCoord ex_range(range_.first, range_.second, p);
        bool representative = (float)range.len() >= consensus_.length() * 0.95f;
        
        const RB_AlleleCoord& range2 = o.seed_ranges_[p2];
        const SeedExt& seed2 = o.seeds_[range2.idx];
        Range range2_ = seed2.getExtendedRange(o.consensus().length());
        RB_AlleleCoord ex_range2(range2_.first, range2_.second, p2);
        bool representative2 = (float)range2.len() >= o.consensus_.length() * 0.95f;
        
        seed_i = range.idx;
        seed_j = range2.idx;
        
        const size_t relax = 5;
        if(representative && representative2) {
            if(ex_range.overlap_len(ex_range2) > 0) {
                if(ex_range.contain(ex_range2, relax)) {
                    contain = true;
                    left = true;
                } else if(ex_range2.contain(ex_range, relax)) {
                    contain = true;
                    left = false;
                } else {
                    left = ex_range.left <= ex_range2.left;
                }
                return true;
            }
        } else if(representative) {
            if(range2.contain(ex_range, relax)) {
                contain = true;
                left = false;
                return true;
            }
        } else if(representative2) {
            if(range.contain(ex_range2, relax)) {
                contain = true;
                left = true;
                return true;
            }
        }

        if(range.right <= range2.right) p++;
        if(range2.right <= range.right) p2++;
    }
    return false;
}

void RB_Repeat::showInfo(const RepeatParameter& rp,
                         CoordHelper& coordHelper) const
{
    cerr << "\trepeat id: " << repeat_id_ << endl;
    cerr << "\tnumber of alleles: " << seeds_.size() << endl;
    cerr << "\tconsensus length: " << consensus_.length() << endl;
    cerr << consensus_ << endl;
    for(size_t i = 0; i < seeds_.size(); i++) {
        const SeedExt& seed = seeds_[i];
        size_t ext_len = seed.getLength();
        if(ext_len < rp.min_repeat_len) continue;
        bool sense_strand = seed.pos.first < coordHelper.forward_length();
        
        cerr << "\t\t" << setw(4) << i << " " << setw(4) << ext_len;
        cerr << " " << (sense_strand ? '+' : '-');
        cerr << " " << setw(10) << seed.pos.first << "  " << setw(10) << seed.pos.second;
        cerr << " " << setw(4) << seed.consensus_pos.first << "  " << setw(4) << seed.consensus_pos.second;
        cerr << " " << (seed.aligned ? "aligned" : "unaligned");
        cerr << endl;
    }
}

template<typename TStr>
void RB_Repeat::generateSNPs(const RepeatParameter& rp, const TStr& s, TIndexOffU grp_id) {
    EList<SeedExt>& seeds = this->seeds();
    
    for(size_t i = 0; i < seeds.size(); i++) {
        SeedExt& seed = seeds[i];
        
        if(!seed.aligned) {
            continue;
        }
        
        if(seed.getLength() < rp.min_repeat_len) {
            continue;
        }
        assert_eq(seed.getLength(), seed.pos.second - seed.pos.first);
        seed.generateSNPs(s, consensus(), snps());
    }
    if(snps().size() > 0) {
        sort(snps_.begin(), snps().end(), SeedSNP::cmpSeedSNPByPos);
    }
}

void RB_Repeat::saveSNPs(ofstream &fp, TIndexOffU grp_id, TIndexOffU& snp_id_base)
{
    string chr_name = "rep" + to_string(repeat_id_);
    
    for(size_t j = 0; j < snps_.size(); j++) {
        SeedSNP& snp = *snps_[j];
        // assign SNPid
        snp.id = (snp_id_base++);
        
        fp << "rps" << snp.id;
        fp << "\t";
        
        if(snp.type == EDIT_TYPE_MM) {
            fp << "single";
        } else if(snp.type == EDIT_TYPE_READ_GAP) {
            fp << "deletion";
        } else if(snp.type == EDIT_TYPE_REF_GAP) {
            fp << "insertion";
        } else {
            assert(false);
        }
        fp << "\t";
        
        fp << chr_name;
        fp << "\t";
        
        fp << snp.pos;
        fp << "\t";
        if(snp.type == EDIT_TYPE_MM || snp.type == EDIT_TYPE_REF_GAP) {
            fp << snp.base;
        } else if(snp.type == EDIT_TYPE_READ_GAP) {
            fp << snp.len;
        } else {
            assert(false);
        }
        fp << endl;
    }
}

float RB_RepeatExt::mergeable(const RB_Repeat& o) const
{
    const EList<RB_AlleleCoord>& ranges = seed_ranges();
    const EList<RB_AlleleCoord>& ranges2 = o.seed_ranges();
    size_t num_overlap_bp = 0;
    size_t p = 0, p2 = 0;
    while(p < ranges.size() && p2 < ranges2.size()) {
        const RB_AlleleCoord& range = ranges[p];
        const RB_AlleleCoord& range2 = ranges2[p2];
        TIndexOffU overlap = range.overlap_len(range2);
        num_overlap_bp += overlap;
        if(range.right <= range2.right) p++;
        else                            p2++;
    }
    size_t num_total_bp = 0, num_total_bp2 = 0;
    for(size_t r = 0; r < ranges.size(); r++) num_total_bp += (ranges[r].right - ranges[r].left);
    for(size_t r = 0; r < ranges2.size(); r++) num_total_bp2 += (ranges2[r].right - ranges2[r].left);
    float portion = float(num_overlap_bp) / float(min(num_total_bp, num_total_bp2));
    return portion;
}

inline void get_next_range(const EList<int>& offsets, size_t i, Range& r, float& avg)
{
    r.first = i;
    for(; r.first < offsets.size() && offsets[r.first] < 0; r.first++);
    r.second = r.first + 1;
    avg = 0.0f;
    if(r.first < offsets.size()) {
        float diff = (float)offsets[r.first] - (float)r.first;
        avg += diff;
    }
    for(; r.second < offsets.size() &&
        offsets[r.second] >= 0 &&
        (r.second == 0 || offsets[r.second] >= offsets[r.second - 1]);
        r.second++) {
        float diff = (float)offsets[r.second] - (float)r.second;
        avg += diff;
    }
    avg /= (float)(r.second - r.first);
    return;
}

template<typename TStr>
bool RB_RepeatExt::align(const RepeatParameter& rp,
                         const TStr& ref,
                         const string& s,
                         const EList<pair<size_t, size_t> >& s_kmer_table,
                         const string& s2,
                         EList<int>& offsets,
                         size_t k,
                         SeedExt& seed,
                         int consensus_approx_left,
                         int consensus_approx_right,
                         size_t left,
                         size_t right,
                         bool debug)
{
    RB_Repeat::seed_merge_tried++;
    
    seed.reset();
    seed.pos.first = left;
    seed.pos.second = right;
    seed.orig_pos.first = seed.pos.first;
    seed.orig_pos.second = seed.orig_pos.first;
    seed.consensus_pos.first = 0;
    seed.consensus_pos.second = consensus_.length();
    seed.aligned = false;
    
    {
        const int query_len = right - left;
        const int consensus_len = consensus_approx_right - consensus_approx_left;
        const int abs_gap_len = max<int>(abs(consensus_len - query_len), 5);
        
        offsets.resize(s2.length());
        offsets.fill(-1);
        pair<size_t, size_t> kmer(0, 0);
        for(size_t i = 0; i + k <= s2.length(); i++) {
            if(i == 0) {
                kmer.first = extract_kmer(s2, i, k);
            } else {
                kmer.first = next_kmer(kmer.first, s2[i+k-1], k);
            }
            size_t lb = s_kmer_table.bsearchLoBound(kmer);
            while(lb < s_kmer_table.size() && kmer.first == s_kmer_table[lb].first) {
                int expected = (int)i + consensus_approx_left;
                int real = (int)s_kmer_table[lb].second;
                int abs_diff = abs(expected - real);
                if(abs_diff <= abs_gap_len * 2 || debug) {
                    if(offsets[i] == -1) {
                        offsets[i] = (int)s_kmer_table[lb].second;
                    } else if(offsets[i] >= 0) {
                        offsets[i] = -2;
                    } else {
                        assert_lt(offsets[i], -1);
                        offsets[i] -= 1;
                    }
                }
                lb++;
            }
            if(offsets[i] > 0 && i + k == s2.length()) {
                for(size_t j = i + 1; j < s2.length(); j++) {
                    offsets[j] = offsets[j-1] + 1;
                }
            }
        }
    }
    
    if(debug) {
        cerr << "initial offsets" << endl;
        for(size_t j = 0; j < offsets.size(); j++) {
            cout << j << ": " << offsets[j] << " " << s2[j] << ": " << (offsets[j] < 0 ? ' ' : s[offsets[j]]) << endl;
        }
    }
    
    // remove inconsistent positions
    Range range; float range_avg;
    get_next_range(offsets, 0, range, range_avg);
    while(range.second < offsets.size()) {
        Range range2; float range_avg2;
        get_next_range(offsets, range.second, range2, range_avg2);
        if(range2.first >= offsets.size())
            break;
        
        assert_leq(range.second, range2.first);
        
        float abs_diff = abs(range_avg - range_avg2);
        if(offsets[range.second - 1] > offsets[range2.first] ||
           (abs_diff > 10.0f && (range.second - range.first < 5 || range2.second - range2.first < 5)) ||
           abs_diff > 20.0f) {
            if(range.second - range.first < range2.second - range2.first) {
                for(size_t i = range.first; i < range.second; i++) {
                    offsets[i] = -1;
                }
                range = range2;
                range_avg = range_avg2;
            } else {
                for(size_t i = range2.first; i < range2.second; i++) {
                    offsets[i] = -1;
                }
            }
        } else {
            range = range2;
            range_avg = range_avg2;
        }
    }
    
    bool weighted_avg_inited = false;
    float weighted_avg = -1.0f;
    for(size_t i = 0; i < offsets.size(); i++) {
        if(offsets[i] < 0)
            continue;

        float diff = (float)offsets[i] - (float)i;
        if(weighted_avg_inited) {
            if(abs(diff - weighted_avg) > 20.0f) {
                offsets[i] = -1;
                continue;
            }
        }
        
        if(weighted_avg < 0.0f) {
            weighted_avg = diff;
            weighted_avg_inited = true;
        } else {
            weighted_avg = 0.8f * weighted_avg + 0.2f * diff;
        }
    }
    
    if(debug) {
        cerr << "after filtering" << endl;
        for(size_t j = 0; j < offsets.size(); j++) {
            cout << j << ": " << offsets[j] << " " << s2[j] << ": " << (offsets[j] < 0 ? ' ' : s[offsets[j]]) << endl;
        }
    }
    
    int i = 0;
    while(i < offsets.size()) {
        // skip non-negative offsets
        for(; i < offsets.size() && offsets[i] >= 0; i++);
        if(i >= offsets.size()) break;
        int j = i;
        for(; j < offsets.size(); j++) {
            if(offsets[j] >= 0)
                break;
        }
        assert(i >= offsets.size() || offsets[i] < 0);
        assert(j >= offsets.size() || offsets[j] >= 0);
        if(i > 0 && j < offsets.size()) {
            i -= 1;
            int left = offsets[i], right = offsets[j];
            assert_geq(left, 0);
            if(left > right) return false;
            assert_leq(left, right);
            int ref_len = right - left + 1;
            int query_len = j - i + 1;
            if(query_len == ref_len) { // match
                for(size_t i2 = i + 1; i2 < j; i2++)
                    offsets[i2] = offsets[i2-1] + 1;
            } else { // deletion or insertion
                bool del = query_len < ref_len;
                size_t gap_len = del ? ref_len - query_len : query_len - ref_len;
                size_t max_len = max(ref_len, query_len);
                const size_t max_mm = max_len / 25 + 1;
                const size_t very_max_mm = max<size_t>(max_len / 2, max_mm);
                ca_s_ = s.substr(left, max_len);
                ca_s2_ = s2.substr(i, max_len);
                
                size_t gap_pos = 0, mm = very_max_mm + 1;
                find_gap_pos(ca_s_,
                             ca_s2_,
                             ca_ed_,
                             ca_ed2_,
                             del,
                             gap_len,
                             very_max_mm,
                             gap_pos,
                             mm,
                             debug);
                
                if(mm > very_max_mm) {
                    return false;
                }
                
                assert_lt(gap_pos, query_len);
                if(del) {
                    for(size_t i2 = i + 1; i2 < j; i2++) {
                        if(i2 - i == gap_pos) {
                            offsets[i2] = offsets[i2-1] + gap_len;
                        } else {
                            offsets[i2] = offsets[i2-1] + 1;
                        }
                    }
                } else {
                    for(size_t i2 = i + 1; i2 < j; i2++) {
                        if(i2 - i >= gap_pos && i2 - i < gap_pos + gap_len) {
                            offsets[i2] = offsets[i2-1];
                        } else {
                            offsets[i2] = offsets[i2-1] + 1;
                        }
                    }
                }
            }
        }
        
        i = j;
    }
    
    if(debug) {
        cerr << "final offsets" << endl;
        for(size_t j = 0; j < offsets.size(); j++) {
            cout << j << ": " << offsets[j] << " " << s2[j] << ": " << (offsets[j] < 0 ? ' ' : s[offsets[j]]) << endl;
        }
    }
    
#ifndef NDEBUG
    for(size_t i = 1; i < offsets.size(); i++) {
        if(offsets[i-1] < 0 || offsets[i] < 0) continue;
        assert_leq(offsets[i-1], offsets[i]);
    }
#endif
    
    size_t b, e;
    for(b = 0; b < offsets.size() && offsets[b] < 0; b++);
    if(b >= offsets.size()) return false;
    assert_lt((size_t)b, offsets.size());
    for(e = offsets.size() - 1; e > b && offsets[e] < 0; e--);
    if(b == e) return false;
    for(size_t p = b; p <= e; p++) {
        if(offsets[p] < 0) return false;
    }

    // try to fill the ends
    if(b > 0) {
        int mm = 0, pb = (int)b;
        for(int i = pb - 1; i >= 0; i--) {
            assert_geq(offsets[i+1], 0);
            if(offsets[i+1] == 0) break;
            if(s2[i] != s[offsets[i+1] - 1])
                mm++;
            if(pb - i < 25 * (mm - 1))
                break;
            offsets[i] = offsets[i+1] - 1;
            b = (size_t)i;
        }
    }
    if(e + 1 < offsets.size()) {
        int mm = 0, prev_end = (int)e;
        for(int i = prev_end + 1; i < offsets.size(); i++) {
            assert_geq(offsets[i-1], 0);
            if(offsets[i-1] + 1 >= s.length()) break;
            if(s2[i] != s[offsets[i-1] + 1])
                mm++;
            if(i - prev_end < 25 * (mm - 1))
                break;
            offsets[i] = offsets[i-1] + 1;
            e = (size_t)i;
        }
    }
    
    assert_geq(seed.pos.first, (TIndexOffU)b);
    seed.pos.first += (TIndexOffU)b;
    seed.pos.second = seed.pos.first + e - b + 1;
    seed.orig_pos.first = seed.pos.first;
    seed.orig_pos.second = seed.pos.first;
    seed.consensus_pos.first = offsets[b];
    seed.consensus_pos.second = offsets[e] + 1;
    seed.aligned = true;
    
    for(int p = b; p < e; p++) {
        assert_geq(offsets[p], 0);
        assert_leq(offsets[p], offsets[p+1]);
        if(offsets[p] + 1 == offsets[p+1]) { // match
            continue;
        } else if(offsets[p] + 1 < offsets[p+1]) { // deletion
            seed.right_gaps.expand();
            seed.right_gaps.back().first = p + 1 - b;
            seed.right_gaps.back().second = offsets[p+1] - offsets[p] - 1;
        } else {
            assert_eq(offsets[p], offsets[p+1]);
            int p2 = p + 1;
            for(; offsets[p2] == offsets[p2+1]; p2++);
            seed.right_gaps.expand();
            seed.right_gaps.back().first = p + 1 - b;
            seed.right_gaps.back().second = (int)p - (int)p2;
            p = p2;
        }
    }
    
    if(debug) {
        string consensus_str = consensus_.substr(seed.consensus_pos.first, seed.consensus_pos.second - seed.consensus_pos.first);
        string seed_str;
        seed.getExtendedSeedSequence(ref, seed_str);
        assert_eq(consensus_str.length(), seed_str.length());
        cerr << "consensus: " << consensus_str << endl;
        cerr << "seed     : " << seed_str << endl;
    }
    
    RB_Repeat::seed_merged++;
    
    return true;
}

bool RB_RepeatExt::isSelfRepeat(const RepeatParameter& rp,
                                const string& s,
                                const EList<pair<size_t, size_t> >& s_kmer_table,
                                EList<int>& offsets,
                                size_t k,
                                bool debug)
{
    offsets.resize(s.length());
    offsets.fill(-1);
    pair<size_t, size_t> kmer(0, 0);
    for(size_t i = 0; i + k <= s.length(); i++) {
        if(i == 0) {
            kmer.first = extract_kmer(s, i, k);
        } else {
            kmer.first = next_kmer(kmer.first, s[i+k-1], k);
        }
        size_t lb = s_kmer_table.bsearchLoBound(kmer);
        while(lb < s_kmer_table.size() && kmer.first == s_kmer_table[lb].first) {
            if(offsets[i] == -1) {
                offsets[i] = (int)s_kmer_table[lb].second;
            } else if(offsets[i] >= 0) {
                offsets[i] = -2;
            } else {
                assert_lt(offsets[i], -1);
                offsets[i] -= 1;
            }
            lb++;
        }
        if(offsets[i] > 0 && i + k == s.length()) {
            for(size_t j = i + 1; j < s.length(); j++) {
                offsets[j] = offsets[j-1] + 1;
            }
        }
    }
    
    if(debug) {
        cerr << "offsets" << endl;
        for(size_t j = 0; j < offsets.size(); j++) {
            cout << j << ": " << offsets[j] << " " << s[j] << ": " << (offsets[j] < 0 ? ' ' : s[offsets[j]]) << endl;
        }
    }
    
    const size_t min_repeat_size = 100;
    size_t repeat_count = 0;
    for(size_t i = 0; i + min_repeat_size < offsets.size(); i++) {
        if(offsets[i] >= -1)
            continue;
        size_t j = i + 1;
        for(; j < offsets.size() && offsets[j] <= -2; j++);
        if(j - i >= min_repeat_size)
            repeat_count++;
        i = j;
    }

    return repeat_count >= 2;
}

template<typename TStr>
bool RB_RepeatExt::merge(const RepeatParameter& rp,
                         const TStr& s,
                         RB_SWAligner& swalginer,
                         const RB_RepeatExt& o,
                         bool contain,
                         size_t seed_i,
                         size_t seed_j,
                         bool debug)
{
    // construct a new consensus sequence
    string prev_consensus = consensus_;
    size_t consensus_add_len = 0;
    {
        assert_lt(seed_i, seeds_.size());
        assert_lt(seed_j, o.seeds_.size());
        const SeedExt& seed = seeds_[seed_i];
        const SeedExt& oseed = o.seeds_[seed_j];
        
        Range range = seed.getExtendedRange(consensus_.length());
        Range orange = oseed.getExtendedRange(o.consensus_.length());
        assert_leq(range.first, orange.first + 10);
        
        consensus_add_len = (int)orange.first - (int)range.first;
        if(!contain) {
            if(range.second <= orange.first) {
                cerr << "something wrong: " << range.first << "-" << range.second << " ";
                cerr << orange.first << "-" << orange.second << " " << o.consensus_.length() << endl;
                assert(false);
            }
            if(range.second > orange.first &&
               range.second - orange.first < o.consensus_.length()) {
                consensus_ += o.consensus_.substr(range.second - orange.first);
            }
        }
    }
    
    EList<pair<size_t, size_t> > merge_list;
    merge_list.reserveExact(o.seed_ranges_.size());
    size_t p = 0, p2 = 0;
    const size_t relax = 5;
    while(p < seed_ranges_.size() && p2 < o.seed_ranges_.size()) {
        const RB_AlleleCoord& range = seed_ranges_[p];
        assert_lt(range.idx, seeds_.size());
        const RB_AlleleCoord& range2 = o.seed_ranges_[p2];
        assert_lt(range2.idx, o.seeds_.size());
        if(range.contain(range2, relax)) {
            merge_list.expand();
            merge_list.back().first = p;
            merge_list.back().second = p2;
            if(debug) {
                cerr << p << ":" << range.left << "-" << range.right << " > ";
                cerr << p2 << ":" << range2.left << "-" << range2.right << endl;
            }
        } else if(range2.contain(range, relax)) {
            merge_list.expand();
            merge_list.back().first = p;
            merge_list.back().second = p2;
            if(debug) {
                cerr << p << ":" << range.left << "-" << range.right << " < ";
                cerr << p2 << ":" << range2.left << "-" << range2.right << endl;
            }
        } else {
            TIndexOffU overlap_len = range.overlap_len(range2);
            bool stored = !merge_list.empty() && merge_list.back().first == p;
            bool stored2 = !merge_list.empty() && merge_list.back().second == p2;
            if(overlap_len > 0) {
                if(!stored && !stored2) {
                    merge_list.expand();
                    merge_list.back().first = p;
                    merge_list.back().second = p2;
                    
                    if(debug) {
                        cerr << p << ":" << range.left << "-" << range.right << " ol ";
                        cerr << p2 << ":" << range2.left << "-" << range2.right << endl;
                    }
                }
            } else {
                if(range2.right <= range.left) {
                    if(!stored2) {
                        merge_list.expand();
                        merge_list.back().first = numeric_limits<size_t>::max();
                        merge_list.back().second = p2;
                        
                        if(debug) {
                            cerr << p << ":" << range.left << "-" << range.right << " <> ";
                            cerr << p2 << ":" << range2.left << "-" << range2.right << endl;
                        }
                    }
                } else {
                    assert_leq(range.right, range2.left);
                    if(!stored) {
                        merge_list.expand();
                        merge_list.back().first = p;
                        merge_list.back().second = numeric_limits<size_t>::max();
                        
                        if(debug) {
                            cerr << p << ":" << range.left << "-" << range.right << " <> ";
                            cerr << p2 << ":" << range2.left << "-" << range2.right << endl;
                        }
                    }
                }
            }
        }
        if(range.right <= range2.right) p++;
        if(range2.right <= range.right) p2++;
    }
    assert(p == seeds_.size() || p2 == o.seeds_.size());
    
    while(p < seed_ranges_.size()) {
        bool stored = !merge_list.empty() && merge_list.back().first == p;
        if(!stored) {
            const RB_AlleleCoord& range = seed_ranges_[p];
            merge_list.expand();
            merge_list.back().first = p;
            merge_list.back().second = numeric_limits<size_t>::max();
            
            if(debug) {
                cerr << p << ":" << range.left << "-" << range.right << " : <> " << endl;
            }
        }
        
        p++;
    }
    
    while(p2 < o.seed_ranges_.size()) {
        bool stored2 = !merge_list.empty() && merge_list.back().second == p2;
        if(!stored2) {
            const RB_AlleleCoord& range2 = o.seed_ranges_[p2];
            merge_list.expand();
            merge_list.back().first = numeric_limits<size_t>::max();
            merge_list.back().second = p2;
            
            if(debug) {
                cerr << ": <> " << p2 << ":" << range2.left << "-" << range2.right << endl;
            }
        }
        
        p2++;
    }
    
    assert(!merge_list.empty());
    
    if(debug) {
        for(size_t i = 0; i < merge_list.size(); i++) {
            cerr << "merge list:" << endl;
            cerr << "\t" << (merge_list[i].first < seed_ranges_.size() ? (int)merge_list[i].first : -1);
            cerr << " " << (merge_list[i].second < o.seed_ranges_.size() ? (int)merge_list[i].second : -1) << endl;
        }
    }
    
    const size_t kmer_len = 12;
    EList<pair<size_t, size_t> > kmer_table;
    build_kmer_table(consensus_, kmer_table, kmer_len);
    EList<int> offsets;
    
    bool self_repeat = isSelfRepeat(rp,
                                    consensus_,
                                    kmer_table,
                                    offsets,
                                    kmer_len,
                                    debug);
    if(self_repeat) {
        consensus_ = prev_consensus;
        return false;
    }
    
    string seq;
    for(size_t i = 0; i < merge_list.size(); i++) {
        size_t seed_id = (merge_list[i].first < seed_ranges_.size() ? seed_ranges_[merge_list[i].first].idx : merge_list[i].first);
        size_t oseed_id = (merge_list[i].second < o.seed_ranges_.size() ? o.seed_ranges_[merge_list[i].second].idx : merge_list[i].second);
        
        if(seed_id < seeds_.size()) {
            if(oseed_id >= o.seeds_.size())
                continue;
            else if(seed_ranges_[merge_list[i].first].contain(o.seed_ranges_[merge_list[i].second]))
                continue;
        }
        
        const SeedExt* seed = (seed_id < seeds_.size() ? &seeds_[seed_id] : NULL);
        const SeedExt* oseed = (oseed_id < o.seeds_.size() ? &o.seeds_[oseed_id] : NULL);
        assert(seed != NULL || oseed != NULL);
        
        size_t left = (seed != NULL ? seed->pos.first : numeric_limits<size_t>::max());
        size_t right = (seed != NULL ? seed->pos.second : 0);
        int consensus_approx_left = (seed != NULL ? seed->consensus_pos.first : numeric_limits<int>::max());
        int consensus_approx_right = (seed != NULL ? seed->consensus_pos.second : 0);
        if(oseed != NULL) {
            if(oseed->pos.first < left) {
                left = oseed->pos.first;
            }
            if(oseed->pos.second > right) {
                right = oseed->pos.second;
            }
            
            int oconsensus_approx_left = oseed->consensus_pos.first + consensus_add_len;
            if(oconsensus_approx_left < consensus_approx_left) {
                consensus_approx_left = oconsensus_approx_left;
            }
            int oconsensus_approx_right = oseed->consensus_pos.second + consensus_add_len;
            if(oconsensus_approx_right > consensus_approx_right) {
                consensus_approx_right = oconsensus_approx_right;
            }
        }

        getString(s, left, right - left, seq);
        if(seed_id >= seeds_.size()) {
            seed_id = seeds_.size();
            seeds_.expand();
        }
        
        /* bool succ = */ align(rp,
                                s,
                                consensus_,
                                kmer_table,
                                seq,
                                offsets,
                                kmer_len,
                                seeds_[seed_id],
                                consensus_approx_left,
                                consensus_approx_right,
                                left,
                                right,
                                debug && right - left == 1080);
    }
    
    while(true) {
        internal_update();
        
        size_t remove_count = 0;
        for(size_t i = 0; i + 1 < seed_ranges_.size(); i++) {
            RB_AlleleCoord& range = seed_ranges_[i];
            if(range.left == numeric_limits<size_t>::max())
                break;
            for(size_t j = i + 1; j < seed_ranges_.size(); j++) {
                RB_AlleleCoord& range2 = seed_ranges_[j];
                if(range2.left == numeric_limits<size_t>::max())
                    break;
                if(range.right <= range2.left)
                    break;
                
                getString(s, range.left, range2.right - range.left, seq);
                
                /* bool succ = */ align(rp,
                                        s,
                                        consensus_,
                                        kmer_table,
                                        seq,
                                        offsets,
                                        kmer_len,
                                        seeds_[range.idx],
                                        seeds_[range.idx].consensus_pos.first,
                                        seeds_[range2.idx].consensus_pos.second,
                                        range.left,
                                        range2.right,
                                        debug && range.left == 692422);
                
                range.left = seeds_[range.idx].pos.first;
                range.right = seeds_[range.idx].pos.second;
                seeds_[range2.idx].reset();
                range2.left = numeric_limits<size_t>::max();
                remove_count++;
            }
        }
        
        if(remove_count <= 0) break;
        
        sort(seeds_.begin(), seeds_.end(), seedCmp);
        seeds_.resize(seeds_.size() - remove_count);
    }
    
    return true;
}

#define DMAX std::numeric_limits<double>::max()

RB_SWAligner::RB_SWAligner()
{
    rnd_.init(0);
}

RB_SWAligner::~RB_SWAligner()
{
    if(sc_) {
        delete sc_;
    }
}

void RB_SWAligner::init_dyn(const RepeatParameter& rp)
{
    const int MM_PEN = 3;
    // const int MM_PEN = 6;
    const int GAP_PEN_LIN = 2;
    // const int GAP_PEN_LIN = (((MM_PEN) * rpt_edit_ + 1) * 1.0);
    const int GAP_PEN_CON = 4;
    // const int GAP_PEN_CON = (((MM_PEN) * rpt_edit_ + 1) * 1.0);
    const int MAX_PEN = MAX_I16;
    
    scoreMin_.init(SIMPLE_FUNC_LINEAR, rp.max_edit * MM_PEN * -1.0, 0.0);
    nCeil_.init(SIMPLE_FUNC_LINEAR, 0.0, 0.0);
    
    penCanIntronLen_.init(SIMPLE_FUNC_LOG, -8, 1);
    penNoncanIntronLen_.init(SIMPLE_FUNC_LOG, -8, 1);
    
    sc_ = new Scoring(
                      DEFAULT_MATCH_BONUS,     // constant reward for match
                      DEFAULT_MM_PENALTY_TYPE,     // how to penalize mismatches
                      MM_PEN,      // max mm penalty
                      MM_PEN,      // min mm penalty
                      MAX_PEN,      // max sc penalty
                      MAX_PEN,      // min sc penalty
                      scoreMin_,       // min score as function of read len
                      nCeil_,          // max # Ns as function of read len
                      DEFAULT_N_PENALTY_TYPE,      // how to penalize Ns in the read
                      DEFAULT_N_PENALTY,           // constant if N pelanty is a constant
                      DEFAULT_N_CAT_PAIR,    // whether to concat mates before N filtering
                      
                      
                      GAP_PEN_CON, // constant coeff for read gap cost
                      GAP_PEN_CON, // constant coeff for ref gap cost
                      GAP_PEN_LIN, // linear coeff for read gap cost
                      GAP_PEN_LIN, // linear coeff for ref gap cost
                      1 /* gGapBarrier */    // # rows at top/bot only entered diagonally
                      );
}

void RB_SWAligner::makePadString(const string& ref,
                                 const string& read,
                                 string& pad,
                                 size_t len)
{
    pad.resize(len);
    
    for(size_t i = 0; i < len; i++) {
        // shift A->C, C->G, G->T, T->A
        pad[i] = "CGTA"[asc2dna[ref[i]]];
        
        if(read[i] == pad[i]) {
            // shift
            pad[i] = "CGTA"[asc2dna[pad[i]]];
        }
    }
    
    int head_len = len / 2;
    size_t pad_start = len - head_len;
    
    for(size_t i = 0; i < head_len; i++) {
        if(read[i] == pad[pad_start + i]) {
            // shift
            pad[pad_start + i] = "CGTA"[asc2dna[pad[pad_start + i]]];
        }
    }
}

int RB_SWAligner::alignStrings(const string &ref,
                               const string &read,
                               EList<Edit>& edits,
                               Coord& coord)
{
    // Prepare Strings
    
    // Read -> BTDnaString
    // Ref -> bit-encoded string
    
    //SwAligner swa;
    
    BTDnaString btread;
    BTString btqual;
    BTString btref;
    BTString btref2;
    
    BTDnaString btreadrc;
    BTString btqualrc;
    
    
    string qual = "";
    for(size_t i = 0; i < read.length(); i++) {
        qual.push_back('I');
    }
    
#if 0
    cerr << "REF : " << ref << endl;
    cerr << "READ: " << read << endl;
    cerr << "QUAL: " << qual << endl;
#endif
    
    btread.install(read.c_str(), true);
    btreadrc = btread;
    btreadrc.reverseComp();
    
    btqual.install(qual.c_str());
    btqualrc = btqual;
    
    btref.install(ref.c_str());
    
    TAlScore min_score = sc_->scoreMin.f<TAlScore >((double)btread.length());
    
    btref2 = btref;
    
    size_t nceil = 0;
    // size_t nrow = btread.length();
    
    // Convert reference string to mask
    for(size_t i = 0; i < btref2.length(); i++) {
        if(toupper(btref2[i]) == 'N') {
            btref2.set(16, i);
        } else {
            int num = 0;
            int alts[] = {4, 4, 4, 4};
            decodeNuc(toupper(btref2[i]), num, alts);
            assert_leq(num, 4);
            assert_gt(num, 0);
            btref2.set(0, i);
            for(int j = 0; j < num; j++) {
                btref2.set(btref2[i] | (1 << alts[j]), i);
            }
        }
    }
    
    
    bool fw = true;
    uint32_t refidx = 0;
    
    swa.initRead(
                 btread,     // read sequence
                 btreadrc,
                 btqual,     // read qualities
                 btqualrc,
                 0,          // offset of first character within 'read' to consider
                 btread.length(), // offset of last char (exclusive) in 'read' to consider
                 *sc_);      // local-alignment score floor
    
    DynProgFramer dpframe(false);
    size_t readgaps = 0;
    size_t refgaps = 0;
    size_t maxhalf = 0;
    
    DPRect rect;
    dpframe.frameSeedExtensionRect(
                                   0,              // ref offset implied by seed hit assuming no gaps
                                   btread.length(),    // length of read sequence used in DP table
                                   btref2.length(),    // length of reference
                                   readgaps,   // max # of read gaps permitted in opp mate alignment
                                   refgaps,    // max # of ref gaps permitted in opp mate alignment
                                   (size_t)nceil,  // # Ns permitted
                                   maxhalf, // max width in either direction
                                   rect);      // DP Rectangle
    
    assert(rect.repOk());
    
    size_t cminlen = 2000, cpow2 = 4;
    
    swa.initRef(
                fw,                 // whether to align forward or revcomp read
                refidx,             // reference ID
                rect,               // DP rectangle
                btref2.wbuf(),      // reference strings
                0,                  // offset of first reference char to align to
                btref2.length(),    // offset of last reference char to align to
                btref2.length(),    // length of reference sequence
                *sc_,               // scoring scheme
                min_score,          // minimum score
                true,               // use 8-bit SSE if positions
                cminlen,            // minimum length for using checkpointing scheme
                cpow2,              // interval b/t checkpointed diags; 1 << this
                false,              // triangular mini-fills?
                false               // is this a seed extension?
                );
    
    
    TAlScore best = std::numeric_limits<TAlScore>::min();
    bool found = swa.align(rnd_, best);
    
#ifdef DEBUGLOG
    cerr << "found: " << found << "\t" << best << "\t" << "minsc: " << min_score << endl;
#endif
    
    if (found) {
#ifdef DEBUGLOG
        //cerr << "CP " << "found: " << found << "\t" << best << "\t" << "minsc: " << min_score << endl;
        cerr << "REF : " << ref << endl;
        cerr << "READ: " << read << endl;
#endif
        
        SwResult res;
        int max_match_len = 0;
        res.reset();
        res.alres.init_raw_edits(&rawEdits_);
        
        found = swa.nextAlignment(res, best, rnd_);
        if (found) {
            edits = res.alres.ned();
            //const TRefOff ref_off = res.alres.refoff();
            //const Coord& coord = res.alres.refcoord();
            coord = res.alres.refcoord();
            //assert_geq(genomeHit._joinedOff + coord.off(), genomeHit.refoff());
            
#ifdef DEBUGLOG
            cerr << "num edits: " << edits.size() << endl;
            cerr << "coord: " << coord.off();
            cerr << ", " << coord.ref();
            cerr << ", " << coord.orient();
            cerr << ", " << coord.joinedOff();
            cerr << endl;
            Edit::print(cerr, edits); cerr << endl;
            Edit::printQAlign(cerr, btread, edits);
#endif
            
            max_match_len = getMaxMatchLen(edits, btread.length());
#ifdef DEBUGLOG
            cerr << "max match length: " << max_match_len << endl;
#endif
        }
#ifdef DEBUGLOG
        cerr << "nextAlignment: " << found << endl;
        cerr << "-------------------------" << endl;
#endif
    }
    
    return 0;
}

void RB_SWAligner::doTest(const RepeatParameter& rp,
                          const string& refstr,
                          const string& readstr)
{
    init_dyn(rp);
    
    doTestCase1(refstr,
                readstr,
                rp.max_edit);
}

void RB_SWAligner::doTestCase1(const string& refstr,
                               const string& readstr,
                               TIndexOffU rpt_edit)
{
    cerr << "doTestCase1----------------" << endl;
    EList<Edit> edits;
    Coord coord;
    
    if (refstr.length() == 0 ||
        readstr.length() == 0) {
        return;
    }
    
    EList<Edit> ed;
    
    string pad;
    makePadString(refstr, readstr, pad, 5);
    
    string ref2 = pad + refstr + pad;
    string read2 = pad + readstr + pad;
    alignStrings(refstr, readstr, edits, coord);
    
    size_t left = pad.length();
    size_t right = left + readstr.length();
    
    edits.reserveExact(ed.size());
    for(size_t i = 0; i < ed.size(); i++) {
        if(ed[i].pos >= left && ed[i].pos <= right) {
            edits.push_back(ed[i]);
            edits.back().pos -= left;
        }
    }
    
    
#if 0
    RepeatGroup rg;
    
    rg.edits = edits;
    rg.coord = coord;
    rg.seq = readstr;
    rg.base_offset = 0;
    
    string chr_name = "rep";
    
    cerr << "REF : " << refstr << endl;
    cerr << "READ: " << readstr << endl;
    size_t snpids = 0;
    rg.buildSNPs(snpids);
    rg.writeSNPs(cerr, chr_name); cerr << endl;
#endif
}


template<typename TStr>
RepeatBuilder<TStr>::RepeatBuilder(TStr& s,
                                   const EList<RefRecord>& szs,
                                   const EList<string>& ref_names,
                                   bool forward_only,
                                   const string& filename) :
s_(s),
coordHelper_(s.length(), forward_only ? s.length() : s.length() / 2, szs, ref_names),
forward_only_(forward_only),
filename_(filename),
forward_length_(forward_only ? s.length() : s.length() / 2)
{
	cerr << "RepeatBuilder: " << filename_ << endl;
}

template<typename TStr>
RepeatBuilder<TStr>::~RepeatBuilder()
{
    for(map<size_t, RB_Repeat*>::iterator it = repeat_map_.begin(); it != repeat_map_.end(); it++) {
        delete it->second;
    }
    repeat_map_.clear();
}

template<typename TStr>
void RepeatBuilder<TStr>::readSA(const RepeatParameter& rp,
                                 BlockwiseSA<TStr>& sa)
{
    TIndexOffU count = 0;
    subSA_.init(s_.length() + 1, rp.min_repeat_len, rp.repeat_count);
    
    while(count < s_.length() + 1) {
        TIndexOffU saElt = sa.nextSuffix();
        count++;
        
        if(count && (count % 10000000 == 0)) {
            cerr << "SA count " << count << endl;
        }
        
        if(saElt == s_.length()) {
            assert_eq(count, s_.length() + 1);
            break;
        }
        
        subSA_.push_back(s_, coordHelper_, saElt, count == s_.length());
    }
    
    cerr << "subSA size is " << subSA_.size() << endl;
    subSA_.dump();
#if 0
    for(size_t i = 0; i < subSA_.size(); i++) {
        TIndexOffU joinedOff = subSA_[i];
        fp << setw(10) << joinedOff << " " << getString(s_, joinedOff, rp.seed_len) << endl;
    }
#endif
    cerr << "subSA mem Usage: " << subSA_.getMemUsage() << endl << endl;
}

template<typename TStr>
void RepeatBuilder<TStr>::readSA(const RepeatParameter &rp,
                                 const BitPackedArray &sa)
{
    TIndexOffU count = 0;

    subSA_.init(s_.length() + 1, rp.min_repeat_len, rp.repeat_count);

    for(size_t i = 0; i < sa.size(); i++) {
        TIndexOffU saElt = sa[i];
        count++;

        if(count && (count % 10000000 == 0)) {
            cerr << "RB count " << count << endl;
        }

        if(saElt == s_.length()) {
            assert_eq(count, s_.length() + 1);
            break;
        }

        subSA_.push_back(s_, coordHelper_, saElt, count == s_.length());
    }

    cerr << "subSA size: " << endl;
    subSA_.dump();
#if 0
    for(size_t i = 0; i < subSA_.size(); i++) {
        TIndexOffU joinedOff = subSA_[i];
        fp << setw(10) << joinedOff << " " << getString(s_, joinedOff, rp.seed_len) << endl;
    }
#endif
    cerr << "subSA mem Usage: " << subSA_.getMemUsage() << endl << endl;
}

template<typename TStr>
void RepeatBuilder<TStr>::build(const RepeatParameter& rp)
{
    string rpt_len_str;
    rpt_len_str = to_string(rp.min_repeat_len) + "-" + to_string(rp.max_repeat_len);

    string seed_filename = filename_ + ".rep." + rpt_len_str + ".seed";
    ofstream fp(seed_filename.c_str());

    swaligner_.init_dyn(rp);

    RB_RepeatManager* repeat_manager = new RB_RepeatManager;
    
    EList<RB_RepeatBase> repeatBases;
    subSA_.buildRepeatBase(s_,
                           coordHelper_,
                           rp.max_repeat_len,
                           repeatBases);
    
    size_t repeat_id = 0;
    for(size_t i = 0; i < repeatBases.size(); i++) {
        addRepeatGroup(rp,
                       repeat_id,
                       *repeat_manager,
                       repeatBases[i],
                       fp);
    }
    
    {
        // Build and test minimizer-based k-mer table
        const size_t window = RB_Minimizer<TStr>::default_w;
        const size_t k = RB_Minimizer<TStr>::default_k;
        RB_KmerTable kmer_table;
        EList<string> seqs;
        seqs.reserveExact(repeat_map_.size());
        for(map<size_t, RB_Repeat*>::const_iterator it = repeat_map_.begin(); it != repeat_map_.end(); it++) {
            const RB_Repeat& repeat = *(it->second);
            assert(repeat.satisfy(rp));
            seqs.expand();
            seqs.back() = repeat.consensus();
        }
        kmer_table.build(seqs,
                         window,
                         k);
        kmer_table.dump(cerr);
        cerr << endl;
        
        string query, rc_query;
        EList<pair<uint64_t, size_t> > minimizers;
        size_t total = 0, num_repeat = 0, correct = 0, false_positive = 0, false_negative = 0;
        for(size_t i = 0; i + rp.min_repeat_len <= forward_length_; i += 1000) {
            if(coordHelper_.getEnd(i) != coordHelper_.getEnd(i + rp.min_repeat_len))
                continue;
            query = getString(s_, i, rp.min_repeat_len);
            rc_query = reverseComplement(query);
            
            TIndexOffU idx = subSA_.find_repeat_idx(s_, query);
            const EList<TIndexOffU>& test_repeat_index = subSA_.getRepeatIndex();
            bool repeat = (idx < test_repeat_index.size());
            bool est_repeat = kmer_table.isRepeat(query,
                                                  rc_query,
                                                  minimizers);
            
            total++;
            if(repeat) num_repeat++;
            if(repeat == est_repeat) {
                correct++;
            } else {
                if(est_repeat) {
                    false_positive++;
                } else {
                    false_negative++;
                    assert(false);
                }
            }
        }
        
        cerr << "total: " << total << endl;
        cerr << "repeat: " << num_repeat << endl;
        cerr << "correct: " << correct << endl;
        cerr << "false positive: " << false_positive << endl;
        cerr << "false negative: " << false_negative << endl;
        cerr << endl;
        
        ELList<RB_Alignment> position2D; EList<RB_Alignment> alignments;
        size_t repeat_total = 0, repeat_aligned = 0;
        const EList<TIndexOffU>& test_repeat_index = subSA_.getRepeatIndex();
        size_t interval = 1;
        if(test_repeat_index.size() >= 1000) {
            interval = test_repeat_index.size() / 1000;
        }
        size_t total_alignments = 0, max_alignments = 0;
        string query2, rc_query2;
        for(size_t i = 0; i < test_repeat_index.size(); i += interval) {
            TIndexOffU saElt_idx = test_repeat_index[i];
            TIndexOffU saElt_idx_end = (i + 1 < test_repeat_index.size() ? test_repeat_index[i+1] : subSA_.size());
            TIndexOffU saElt_size = saElt_idx_end - saElt_idx;
            TIndexOffU saElt = subSA_[saElt_idx];
            query = getString(s_, saElt, rp.min_repeat_len);
            query2 = query;
            
            // introduce three mismatches into a query when the query is at lest 100-bp
            if(rp.min_repeat_len >= 100) {
                const size_t mid_pos1 = (size_t)(rp.min_repeat_len * 0.1);
                if(query2[mid_pos1] == 'A') {
                    query2[mid_pos1] = 'C';
                } else {
                    query2[mid_pos1] = 'A';
                }
                const size_t mid_pos2 = (size_t)(rp.min_repeat_len * 0.5);
                if(query2[mid_pos2] == 'C') {
                    query2[mid_pos2] = 'G';
                } else {
                    query2[mid_pos2] = 'C';
                }
                const size_t mid_pos3 = (size_t)(rp.min_repeat_len * 0.9);
                if(query2[mid_pos3] == 'G') {
                    query2[mid_pos3] = 'T';
                } else {
                    query2[mid_pos3] = 'G';
                }
            }
            
            repeat_total += saElt_size;
            
            size_t found = 0;
            size_t cur_alignments = 0;
            if(kmer_table.isRepeat(query2, minimizers)) {
                kmer_table.findAlignments(query2,
                                          minimizers,
                                          position2D,
                                          alignments);
                total_alignments += (alignments.size() * saElt_size);
                cur_alignments = alignments.size();
                TIndexOffU baseoff = 0;
                for(size_t s = 0; s < seqs.size(); s++) {
                    int spos = seqs[s].find(query);
                    if(spos != string::npos) {
                        for(size_t a = 0; a < alignments.size(); a++) {
                            if(alignments[a].pos == baseoff + spos) {
                                found++;
                            }
                        }
                    }
                    baseoff += seqs[s].length();
                }
                
                assert_leq(found, 1);
            }
        
            rc_query = reverseComplement(query);
            rc_query2 = reverseComplement(query2);
            size_t rc_found = 0;
            if(kmer_table.isRepeat(rc_query2, minimizers)) {
                kmer_table.findAlignments(rc_query2,
                                          minimizers,
                                          position2D,
                                          alignments);
                total_alignments += (alignments.size() * saElt_size);
                cur_alignments += alignments.size();
                if(cur_alignments > max_alignments) {
                    max_alignments = cur_alignments;
                }
                
                TIndexOffU baseoff = 0;
                for(size_t s = 0; s < seqs.size(); s++) {
                    int spos = seqs[s].find(rc_query);
                    if(spos != string::npos) {
                        for(size_t a = 0; a < alignments.size(); a++) {
                            if(alignments[a].pos == baseoff + spos) {
                                rc_found++;
                            }
                        }
                    }
                    baseoff += seqs[s].length();
                }
                assert_leq(rc_found, 1);
            }

            if(found + rc_found == 1) {
                repeat_aligned += saElt_size;
            }
        }

        cerr << "num repeats: " << repeat_total << endl;
        cerr << "repeat aligned using minimizers: " << repeat_aligned << endl;
        cerr << "average number of alignments: " << (float)total_alignments / repeat_total << endl;
        cerr << "max alignment: " << max_alignments << endl;
        cerr << endl;
    }
    
    cerr << "number of repeats is " << repeat_map_.size() << endl;
    
    size_t total_rep_seq_len = 0, total_allele_seq_len = 0;
    size_t i = 0;
    for(map<size_t, RB_Repeat*>::iterator it = repeat_map_.begin(); it != repeat_map_.end(); it++, i++) {
        RB_Repeat& repeat = *(it->second);
        if(!repeat.satisfy(rp))
            continue;
        
        repeat.saveSeedExtension(rp,
                                 s_,
                                 coordHelper_,
                                 i,
                                 fp,
                                 total_rep_seq_len,
                                 total_allele_seq_len);
    }
    
    size_t total_qual_seeds = 0;
    for(map<size_t, RB_Repeat*>::iterator it = repeat_map_.begin(); it != repeat_map_.end(); it++, i++) {
        RB_Repeat& repeat = *(it->second);
        EList<SeedExt>& seeds = repeat.seeds();
        for(size_t i = 0; i < seeds.size(); i++) {
            if(seeds[i].getLength() < rp.min_repeat_len)
                continue;
            total_qual_seeds++;
        }
    }
    
    cerr << "total repeat sequence length: " << total_rep_seq_len << endl;
    cerr << "total allele sequence length: " << total_allele_seq_len << endl;
    cerr << "total number of seeds including those that are position-wise different, but sequence-wise identical: " << total_qual_seeds << endl;
    cerr << endl;
    fp << "total repeat sequence length: " << total_rep_seq_len << endl;
    fp << "total allele sequence length: " << total_allele_seq_len << endl;
    fp << "total number of seeds including those that are position-wise different, but sequence-wise identical: " << total_qual_seeds << endl;
    fp.close();
    
    delete repeat_manager;
    repeat_manager = NULL;
    
    // remove non-qualifying repeats
    //  and update repeat IDs for those remaining
    {
        map<size_t, RB_Repeat*> temp_repeat_map;
        size_t i = 0;
        for(map<size_t, RB_Repeat*>::iterator it = repeat_map_.begin(); it != repeat_map_.end(); it++) {
            RB_Repeat* repeat = it->second;
            if(!repeat->satisfy(rp)) {
                delete repeat;
                continue;
            }
            repeat->repeat_id(i);
            temp_repeat_map[repeat->repeat_id()] = repeat;
            i++;
        }
        repeat_map_ = temp_repeat_map;
    }
    
    const bool sanity_check = true;
    if(sanity_check) {
        EList<pair<size_t, size_t> > kmer_table;
        string query;
        size_t total = 0, match = 0;
        EList<TIndexOffU> positions;
        const EList<TIndexOffU>& test_repeat_index = subSA_.getRepeatIndex();
        size_t interval = 1;
        if(test_repeat_index.size() >= 10000) {
            interval = test_repeat_index.size() / 10000;
        }
        for(size_t i = 0; i < test_repeat_index.size(); i += interval) {
            TIndexOffU saElt_idx = test_repeat_index[i];
            TIndexOffU saElt_idx_end = (i + 1 < test_repeat_index.size() ? test_repeat_index[i+1] : subSA_.size());
            positions.clear();
            for(size_t j = saElt_idx; j < saElt_idx_end; j++) {
                positions.push_back(subSA_[j]);
#ifndef NDEBUG
                if(j > saElt_idx) {
                    TIndexOffU lcp_len = getLCP(s_,
                                                coordHelper_,
                                                positions[0],
                                                positions.back(),
                                                rp.min_repeat_len);
                    assert_eq(lcp_len, rp.min_repeat_len);
                }
                
                TIndexOffU saElt = subSA_[j];
                TIndexOffU start = coordHelper_.getStart(saElt);
                TIndexOffU start2 = coordHelper_.getStart(saElt + rp.min_repeat_len - 1);
                assert_eq(start, start2);
#endif
            }
            
            TIndexOffU saElt = subSA_[saElt_idx];
            size_t true_count = saElt_idx_end - saElt_idx;
            getString(s_, saElt, rp.min_repeat_len, query);
            total++;
            
            size_t count = 0, rc_count = 0;
            string rc_query = reverse_complement(query);
            for(map<size_t, RB_Repeat*>::iterator it = repeat_map_.begin(); it != repeat_map_.end(); it++) {
                RB_Repeat& repeat = *(it->second);
                int pos = repeat.consensus().find(query);
                if(pos != string::npos) {
                    for(size_t s = 0; s < repeat.seeds().size(); s++) {
                        SeedExt& seed = repeat.seeds()[s];
                        string seq;
                        seed.getExtendedSeedSequence(s_, seq);
                        if(seq.find(query) != string::npos)
                            count++;
                    }
                }
                pos = repeat.consensus().find(rc_query);
                if(pos != string::npos) {
                    for(size_t s = 0; s < repeat.seeds().size(); s++) {
                        SeedExt& seed = repeat.seeds()[s];
                        string seq;
                        seed.getExtendedSeedSequence(s_, seq);
                        if(seq.find(rc_query) != string::npos)
                            rc_count++;
                    }
                }
            }
            
            if(count == true_count || rc_count == true_count) {
                match++;
            } else if(total - match <= 10) {
                cerr << "   query: " << query << endl;
                cerr << "rc_query: " << rc_query << endl;
                cerr << "true count: " << true_count << endl;
                cerr << "found count: " << count << endl;
                cerr << "rc found count: " << rc_count << endl;
                cerr << endl;
            }
        }
        
        cerr << "RepeatBuilder: sanity check: " << match << " passed (out of " << total << ")" << endl << endl;
    }
}

template<typename TStr>
void RepeatBuilder<TStr>::writeHaploType(const EList<SeedHP>& haplo_lists,
                                         const EList<SeedExt>& seeds,
                                         TIndexOffU& hapl_id_base,
                                         const string& seq_name,
                                         ostream &fp)
{
    if(haplo_lists.size() == 0) {
        return;
    }

    for(size_t i = 0; i < haplo_lists.size(); i++) {
        const SeedHP &haploType = haplo_lists[i];

        fp << "rpht" << hapl_id_base++;
        fp << "\t" << seq_name;
        fp << "\t" << haploType.range.first;
        fp << "\t" << haploType.range.second;
        fp << "\t";

        assert_gt(haploType.snpIDs.size(), 0);

        for (size_t j = 0; j < haploType.snpIDs.size(); j++) {
            if(j) {
                fp << ",";
            }
            fp << haploType.snpIDs[j];
        }
        fp << endl;
    }
}

template<typename TStr>
void RepeatBuilder<TStr>::writeAllele(TIndexOffU grp_id,
                                      TIndexOffU allele_id,
                                      Range range,
                                      const string& seq_name,
                                      TIndexOffU baseoff,
                                      const EList<SeedExt>& seeds,
                                      ostream &fp)
{
    // >rpt_name*0\trep\trep_pos\trep_len\tpos_count\t0
    // chr_name:pos:direction chr_name:pos:direction
    //
    // >rep1*0	rep	0	100	470	0
    // 22:112123123:+ 22:1232131113:+
    //
    size_t snp_size = seeds[range.first].snps.size();
    size_t pos_size = range.second - range.first;
    fp << ">";
    fp << "rpt_" << grp_id << "*" << allele_id;
    fp << "\t" << seq_name;
    fp << "\t" << baseoff + seeds[range.first].consensus_pos.first;
    fp << "\t" << seeds[range.first].consensus_pos.second - seeds[range.first].consensus_pos.first;
    fp << "\t" << pos_size;
    fp << "\t" << snp_size;

    fp << "\t";
    for(size_t i = 0; i < snp_size; i++) {
        if(i > 0) {fp << ",";}
        fp << "rps" << seeds[range.first].snps[i]->id;
    }
    fp << endl;

    // print positions
    for(size_t i = 0; i < pos_size; i++) {
        if(i > 0 && (i % 10 == 0)) {
            fp << endl;
        }

        if(i % 10 != 0) {
            fp << " ";
        }            
        string chr_name;
        TIndexOffU pos_in_chr;
        TIndexOffU joinedOff = seeds[range.first + i].pos.first;

        bool fw = true;
        if(joinedOff < forward_length_) {
            fw = true;
        } else {
            fw = false;
            joinedOff = s_.length() - joinedOff - seeds[range.first].getLength();
        }

        coordHelper_.getGenomeCoord(joinedOff, chr_name, pos_in_chr);

        char direction = fw ? '+' : '-';
        fp << chr_name << ":" << pos_in_chr << ":" << direction;
    }

    fp << endl;
}

template<typename TStr>
void RepeatBuilder<TStr>::writeSNPs(ostream& fp,
                                    const string& rep_chr_name,
                                    const ESet<Edit>& editset)
{
    for(size_t i = 0; i < editset.size(); i++) {
        const Edit& ed = editset[i];

        if(ed.isMismatch()) {
            fp << "rps" << ed.snpID;
            fp << "\t" << "single";
            fp << "\t" << rep_chr_name;
            fp << "\t" << ed.pos;
            fp << "\t" << ed.qchr;
            fp << endl;
        } else {
            assert(false);
        }
    }
}


template<typename TStr>
void RepeatBuilder<TStr>::saveFile(const RepeatParameter& rp)
{
    saveRepeats(rp);
}

bool RB_RepeatManager::checkRedundant(const RepeatParameter& rp,
                                      const map<size_t, RB_Repeat*>& repeat_map,
                                      const EList<TIndexOffU>& positions,
                                      EList<size_t>& to_remove) const
{
    to_remove.clear();
    bool replace = false;
    for(size_t i = 0; i < positions.size(); i++) {
        TIndexOffU seed_pos = positions[i];
        Range seed_range(seed_pos + rp.seed_len, 0);
        map<Range, EList<size_t> >::const_iterator it = range_to_repeats_.upper_bound(seed_range);
        if(it == range_to_repeats_.end())
            continue;
        
        for(map<Range, EList<size_t> >::const_reverse_iterator rit(it); rit != range_to_repeats_.rend(); rit++) {
            Range repeat_range = rit->first;
            if(repeat_range.first + rp.max_repeat_len <= seed_pos)
                break;
            
            const EList<size_t>& repeat_ids = rit->second;
            assert_gt(repeat_ids.size(), 0);
            for(size_t r = 0; r < repeat_ids.size(); r++) {
                size_t repeat_id = repeat_ids[r];
                size_t idx = to_remove.bsearchLoBound(repeat_id);
                if(idx < to_remove.size() && to_remove[idx] == repeat_id)
                    continue;
                
                bool overlap = seed_pos < repeat_range.second && seed_pos + rp.seed_len > repeat_range.first;
                if(!overlap)
                    continue;
                
                map<size_t, RB_Repeat*>::const_iterator it2 = repeat_map.find(repeat_id);
                assert(it2 != repeat_map.end());
                const RB_Repeat& repeat = *(it2->second);
                const EList<RB_AlleleCoord>& allele_ranges = repeat.seed_ranges();
                size_t num_contain = 0, num_overlap = 0, num_close = 0, num_overlap_bp = 0;
                size_t p = 0, p2 = 0;
                while(p < positions.size() && p2 < allele_ranges.size()) {
                    RB_AlleleCoord range;
                    range.left = positions[p];
                    range.right = positions[p] + rp.seed_len;
                    RB_AlleleCoord range2 = allele_ranges[p2];
                    if(range2.contain(range)) {
                        num_contain++;
                        num_overlap_bp += rp.seed_len;
                    } else {
                        TIndexOffU overlap = range2.overlap_len(range);
                        if(overlap > 0) {
                            num_overlap++;
                            num_overlap_bp += (range2.right - range.left);
                        } else if(range.right + 10 > range2.left && range2.right + 10 > range.left) {
                            num_close++;
                        }
                    }
                    if(range.right <= range2.right) p++;
                    else                            p2++;
                }
                
                // if the number of matches is >= 90% of positions in the smaller group
                if((num_contain + num_overlap) * 10 + num_close * 8 >= min(positions.size(), allele_ranges.size()) * 9) {
                    if(positions.size() <= allele_ranges.size()) {
                        return true;
                    } else {
                        replace = true;
                        to_remove.push_back(repeat_id);
                        to_remove.sort();
                    }
                }
            }
        }
        
        // DK - check this out
        if(replace)
            break;
    }
    return false;
}

void RB_RepeatManager::addRepeat(const RB_Repeat* repeat)
{
    const EList<RB_AlleleCoord>& allele_ranges = repeat->seed_ranges();
    for(size_t i = 0; i < allele_ranges.size(); i++) {
        Range allele_range(allele_ranges[i].left, allele_ranges[i].right);
        addRepeat(allele_range, repeat->repeat_id());
    }
}

void RB_RepeatManager::addRepeat(Range range, size_t repeat_id)
{
    if(range_to_repeats_.find(range) == range_to_repeats_.end()) {
        range_to_repeats_[range] = EList<size_t>();
    }
    EList<size_t>& repeat_ids = range_to_repeats_[range];
    size_t idx = repeat_ids.bsearchLoBound(repeat_id);
    if(idx < repeat_ids.size() && repeat_ids[idx] == repeat_id)
        return;
    repeat_ids.push_back(repeat_id);
    repeat_ids.sort();
}

void RB_RepeatManager::removeRepeat(const RB_Repeat* repeat)
{
    const EList<RB_AlleleCoord>& allele_ranges = repeat->seed_ranges();
    for(size_t p = 0; p < allele_ranges.size(); p++) {
        Range range(allele_ranges[p].left, allele_ranges[p].right);
        removeRepeat(range, repeat->repeat_id());
    }
}

void RB_RepeatManager::removeRepeat(Range range, size_t repeat_id)
{
    EList<size_t>& repeat_ids = range_to_repeats_[range];
    TIndexOffU idx = repeat_ids.bsearchLoBound(repeat_id);
    if(idx < repeat_ids.size() && repeat_id == repeat_ids[idx]) {
        repeat_ids.erase(idx);
        if(repeat_ids.empty())
            range_to_repeats_.erase(range);
    }
}

void RB_RepeatManager::showInfo(const RepeatParameter& rp,
                                CoordHelper& coordHelper,
                                const map<size_t, RB_Repeat*>& repeat_map,
                                size_t num_case) const
{
    size_t count = 0;
    for(map<Range, EList<size_t> >::const_iterator it = range_to_repeats_.begin();
        it != range_to_repeats_.end(); it++) {
        const Range& range = it->first;
        if(range.second - range.first < rp.min_repeat_len) continue;
        map<Range, EList<size_t> >::const_iterator jt = it; jt++;
        for(; jt != range_to_repeats_.end(); jt++) {
            const Range& range2 = jt->first;
            if(range2.second - range2.first < rp.min_repeat_len) continue;
            if(range.second <= range2.first)
                break;
            if(count < num_case) {
                cerr << "range (" << range.first << ", " << range.second << ") vs. range2 (";
                cerr << range2.first << ", " << range2.second << ")" << endl;
                for(size_t i = 0; i < it->second.size(); i++) {
                    cerr << "\t1 " << it->second[i] << endl;
                    const RB_Repeat* repeat = repeat_map.find(it->second[i])->second;
                    repeat->showInfo(rp, coordHelper);
                }
                for(size_t i = 0; i < jt->second.size(); i++) {
                    cerr << "\t2 " << jt->second[i] << endl;
                    const RB_Repeat* repeat = repeat_map.find(jt->second[i])->second;
                    repeat->showInfo(rp, coordHelper);
                }
                cerr << endl << endl;
            }
            count++;
        }
    }
    cerr << "ShowInfo - count: " << count << endl;
}

template<typename TStr>
void RepeatBuilder<TStr>::reassignSeeds(const RepeatParameter& rp,
                                        size_t repeat_bid, // repeat begin id
                                        size_t repeat_eid, // repeat end id
                                        EList<SeedExt>& seeds)
{
    assert_lt(repeat_bid, repeat_eid);
    EList<bool> updated;
    updated.resizeExact(repeat_eid - repeat_bid);
    updated.fillZero();
    
    string seq;
    for(size_t s = 0; s < seeds.size(); s++) {
        SeedExt& seed = seeds[s];
        size_t max_repeat_id = repeat_bid;
        size_t max_ext_len = 0;
        for(size_t i = repeat_bid; i < repeat_eid; i++) {
            map<size_t, RB_Repeat*>::iterator it = repeat_map_.find(i);
            assert(it != repeat_map_.end());
            const RB_Repeat& repeat = *(it->second);
            const string& consensus = repeat.consensus();
            const size_t ext_len = (consensus.length() - rp.seed_len) / 2;
            
            if(seed.bound.first + ext_len > seed.pos.first)
                continue;
            if(seed.pos.second + ext_len > seed.bound.second)
                continue;
            
            getString(s_, seed.pos.first - ext_len, rp.seed_len + ext_len * 2, seq);
            assert_eq(consensus.length(), seq.length());
            size_t tmp_ext_len = 0;
            for(size_t j = 0; j < ext_len; j++, tmp_ext_len++) {
                if(seq[ext_len - j - 1] != consensus[ext_len - j - 1] ||
                   seq[ext_len + rp.seed_len + j] != consensus[ext_len + rp.seed_len + j]) {
                    break;
                }
            }
            
            if(tmp_ext_len > max_ext_len) {
                max_repeat_id = i;
                max_ext_len = tmp_ext_len;
            }
        }
        if(rp.seed_len + max_ext_len * 2 >= rp.min_repeat_len) {
            seed.pos.first -= max_ext_len;
            seed.pos.second += max_ext_len;
            map<size_t, RB_Repeat*>::iterator it = repeat_map_.find(max_repeat_id);
            assert(it != repeat_map_.end());
            RB_Repeat& repeat = *(it->second);
            const string& consensus = repeat.consensus();
            const size_t ext_len = (consensus.length() - rp.seed_len) / 2;
            assert_leq(max_ext_len, ext_len);
            seed.consensus_pos.first = ext_len - max_ext_len;
            seed.consensus_pos.second = consensus.length() - (ext_len - max_ext_len);
            assert_leq(seed.consensus_pos.second - seed.consensus_pos.first, consensus.length());
            repeat.addSeed(seed);
            updated[max_repeat_id - repeat_bid] = true;
        }
    }
    
    for(size_t i = repeat_bid; i < repeat_eid; i++) {
        if(!updated[i - repeat_bid])
            continue;
        map<size_t, RB_Repeat*>::iterator it = repeat_map_.find(i);
        assert(it != repeat_map_.end());
        RB_Repeat& repeat = *(it->second);
        repeat.update();
    }
}

/**
 * TODO
 * @brief 
 *
 * @param rpt_seq
 * @param rpt_range
 */
template<typename TStr>
void RepeatBuilder<TStr>::addRepeatGroup(const RepeatParameter& rp,
                                         size_t& repeat_id,
                                         RB_RepeatManager& repeat_manager,
                                         const RB_RepeatBase& repeatBase,
                                         ostream& fp)
{
    RB_Repeat* repeat = new RB_Repeat;
    repeat->repeat_id(repeat_id);
    repeat->init(rp,
                 s_,
                 coordHelper_,
                 subSA_,
                 repeatBase);
    
    assert(repeat_map_.find(repeat->repeat_id()) == repeat_map_.end());
    repeat_map_[repeat->repeat_id()] = repeat;
    repeat_id++;
}

template<typename TStr>
bool RepeatBuilder<TStr>::checkSequenceMergeable(const string& ref,
                                                 const string& read,
                                                 EList<Edit>& edits,
                                                 Coord& coord,
                                                 TIndexOffU rpt_len,
                                                 TIndexOffU max_edit)
{
    size_t max_matchlen = 0;
    EList<Edit> ed;

    string pad;
    swaligner_.makePadString(ref, read, pad, 5);

    string ref2 = pad + ref;
    string read2 = pad + read;

    swaligner_.alignStrings(ref2, read2, ed, coord);

    // match should start from pad string
    if(coord.off() != 0) {
        return false;
    }

    // no edits on pad string
    if(ed.size() > 0 && ed[0].pos < pad.length()) {
        return false;
    }

    size_t left = pad.length();
    size_t right = left + read.length();

    edits.clear();
    edits.reserveExact(ed.size());
    for(size_t i = 0; i < ed.size(); i++) {
        if(ed[i].pos >= left && ed[i].pos <= right) {
            edits.push_back(ed[i]);
            edits.back().pos -= left;
        }
    }

    max_matchlen = getMaxMatchLen(edits, read.length());

#ifdef DEBUGLOG
    {
        cerr << "After pad removed" << endl;
        BTDnaString btread;
        btread.install(read.c_str(), true);
        Edit::print(cerr, edits); cerr << endl;
        Edit::printQAlign(cerr, btread, edits);
    }
#endif

    return (max_matchlen >= rpt_len);
}


template<typename TStr>
void RepeatBuilder<TStr>::saveRepeats(const RepeatParameter &rp)
{
    ios_base::openmode mode = ios_base::out;
    if(rp.append_result) {
        mode |= ios_base::app;
    } else {
        mode |= ios_base::trunc;
    }
    // Generate SNPs
    size_t i = 0;
    for(map<size_t, RB_Repeat*>::iterator it = repeat_map_.begin(); it != repeat_map_.end(); it++, i++) {
        RB_Repeat& repeat = *(it->second);
        if(!repeat.satisfy(rp))
            continue;

        // for each repeats
        repeat.generateSNPs(rp, s_, i);
    }

    // save snp, consensus sequenuce, info
    string snp_fname = filename_ + ".rep.snp";
    string info_fname = filename_ + ".rep.info";
    string hapl_fname = filename_ + ".rep.haplotype";

    ofstream snp_fp(snp_fname.c_str(), mode);
    ofstream info_fp(info_fname.c_str(), mode);
    ofstream hapl_fp(hapl_fname.c_str(), mode);

    const string repName = "rep" + to_string(rp.min_repeat_len) + "-" + to_string(rp.max_repeat_len);
    
    i = 0;
    TIndexOffU consensus_baseoff = 0;
    TIndexOffU snp_id_base = 0;
    TIndexOffU hapl_id_base = 0;
    for(map<size_t, RB_Repeat*>::iterator it = repeat_map_.begin(); it != repeat_map_.end(); it++, i++) {
        RB_Repeat& repeat = *(it->second);
        if(!repeat.satisfy(rp))
            continue;

        // for each repeats
        repeat.saveSNPs(snp_fp, i, snp_id_base);
        saveAlleles(rp,
                    repName,
                    repeat,
                    info_fp,
                    hapl_fp,
                    i,
                    hapl_id_base,
                    consensus_baseoff);

        consensus_baseoff += repeat.consensus().length();
    }

    snp_fp.close();
    info_fp.close();
    hapl_fp.close();

    // save all consensus sequence
    saveConsensus(rp, repName);
}

template<typename TStr>
void RepeatBuilder<TStr>::saveConsensus(const RepeatParameter &rp,
                                        const string& repName) {
    ios_base::openmode mode = ios_base::out;
    if(rp.append_result) {
        mode |= ios_base::app;
    } else {
        mode |= ios_base::trunc;
    }

    string fa_fname = filename_ + ".rep.fa";
    ofstream fa_fp(fa_fname.c_str(), mode);

    fa_fp << ">" << repName << endl;

    size_t oskip = 0;
    for(map<size_t, RB_Repeat*>::iterator it = repeat_map_.begin(); it != repeat_map_.end(); it++) {
        RB_Repeat& repeat = *(it->second);
        if(!repeat.satisfy(rp))
            continue;
        
        // for each repeats
        const string& constr = repeat.consensus();
        size_t si = 0;
        size_t constr_len = constr.length();
        while(si < constr_len) {
            size_t out_len = std::min((size_t)(output_width - oskip), (size_t)(constr_len - si));
            fa_fp << constr.substr(si, out_len);

            if((oskip + out_len) == output_width) {
                fa_fp << endl;
                oskip = 0;
            } else {
                // last line
                oskip = oskip + out_len;
            }

            si += out_len;
        }
    }
    if(oskip) {
        fa_fp << endl;
    }

    fa_fp.close();
}

template<typename TStr>
void RepeatBuilder<TStr>::saveAlleles(
        const RepeatParameter& rp,
        const string& repName,
        RB_Repeat& repeat,
        ofstream& fp,
        ofstream& hapl_fp,
        TIndexOffU grp_id,
        TIndexOffU& hapl_id_base,
        TIndexOffU baseoff)
{
    const EList<SeedExt>& seeds = repeat.seeds();
    EList<SeedHP> haplo_lists;
    Range range(0, seeds.size());

    int allele_id = 0;

    for(size_t sb = range.first; sb < range.second;) {
        size_t se = sb + 1;
        for(; se < range.second; se++) {
            if(!(SeedExt::isSameConsensus(seeds[sb], seeds[se])
                 && SeedExt::isSameSNPs(seeds[sb], seeds[se])
                 && seeds[sb].aligned == seeds[se].aligned)) {
                break;
            }
        }

        if(!seeds[sb].aligned) {
            sb = se;
            continue;
        }

        if(seeds[sb].getLength() < rp.min_repeat_len) {
            sb = se;
            continue;
        }

        // [sb, se) are same alleles
        writeAllele(grp_id, allele_id, Range(sb, se),
                    repName, baseoff, seeds, fp);
        generateHaploType(Range(sb, se), seeds, haplo_lists);

        allele_id++;
        sb = se;
    }

    // sort HaploType List by pos
    haplo_lists.sort();
    writeHaploType(haplo_lists, seeds, hapl_id_base, repName, hapl_fp);
}

template<typename TStr>
void RepeatBuilder<TStr>::generateHaploType(Range range, const EList<SeedExt> &seeds, EList<SeedHP> &haplo_list)
{
    const EList<SeedSNP *>& snps = seeds[range.first].snps;
    if(snps.size() == 0)
        return;

    // right-most position
    TIndexOffU max_right_pos = seeds[range.first].consensus_pos.second - 1;

#ifndef NDEBUG
    for(size_t i = 0; i < snps.size(); i++) {
        assert_leq(snps[i]->pos, max_right_pos);
    }
#endif

    // create haplotypes of at least 16 bp long to prevent combinations of SNPs
    // break a list of SNPs into several haplotypes if a SNP is far from the next SNP in the list
    const TIndexOffU min_ht_len = 16;
    size_t eb = 0, ee = 1;
    while(ee < snps.size() + 1) {
        if(ee == snps.size() ||
           snps[eb]->pos + (min_ht_len << 1) < snps[ee]->pos) {
            TIndexOffU left_pos = snps[eb]->pos;
            TIndexOffU right_pos = snps[ee-1]->pos;
            
            if(snps[ee-1]->type == EDIT_TYPE_READ_GAP) {
                right_pos += snps[ee-1]->len;
            }
            if(left_pos + min_ht_len - 1 > right_pos) {
                right_pos = left_pos + min_ht_len - 1;
            }
            right_pos = min<TIndexOffU>(max_right_pos, right_pos);
            assert_leq(left_pos, right_pos);

            SeedHP seedHP;

            seedHP.range.first = left_pos;
            seedHP.range.second = right_pos;
            for(size_t i = eb; i < ee; i++) {
                string snp_ids = "rps" + to_string(snps[i]->id);
                seedHP.snpIDs.push_back(snp_ids);
            }

            // Add to haplo_list
            bool found = false;
            for(size_t i = 0; i < haplo_list.size(); i++) {
                if(haplo_list[i] == seedHP) {
                    found = true;
                    break;
                }
            }
            if(!found) {
                // Add
                haplo_list.push_back(seedHP);
            }

            eb = ee;
        }
        ee++;
    }
}


void RB_SubSA::init(TIndexOffU sa_size,
                    TIndexOffU seed_len,
                    TIndexOffU seed_count)
{
    assert_gt(sa_size, 1);
    sa_size_ = sa_size;
    assert_gt(seed_len, 0);
    seed_len_ = seed_len;
    assert_gt(seed_count, 1);
    seed_count_ = seed_count;
    temp_suffixes_.clear();

    repeat_list_.clear();
    repeat_index_.clear();
}

template<typename TStr>
void RB_SubSA::push_back(const TStr& s,
                         CoordHelper& coordHelper,
                         TIndexOffU saElt,
                         bool lastInput)
{
    if(saElt + seed_len() <= coordHelper.getEnd(saElt)) {
        if(seed_count_ == 1) {
            repeat_list_.push_back(saElt);
        } else {
            assert_gt(seed_count_, 1);
            
            if(temp_suffixes_.empty()) {
                temp_suffixes_.push_back(saElt);
                return;
            }
            
            TIndexOffU prev_saElt = temp_suffixes_.back();
            
            // calculate common prefix length between two text.
            //   text1 is started from prev_saElt and text2 is started from saElt
            bool same = isSameSequenceUpto(s,
                                           coordHelper,
                                           prev_saElt,
                                           saElt,
                                           seed_len_);
            
            if(same) {
                temp_suffixes_.push_back(saElt);
            }
            
            if(!same || lastInput) {
                if(temp_suffixes_.size() >= seed_count_) {
                    repeat_index_.push_back(repeat_list_.size());
                    temp_suffixes_.sort();
                    for(size_t pi = 0; pi < temp_suffixes_.size(); pi++) {
                        repeat_list_.push_back(temp_suffixes_[pi]);
                    }
                }
                temp_suffixes_.clear();
                if(!lastInput) {
                    temp_suffixes_.push_back(saElt);
                } else {
                    temp_suffixes_.nullify();
                }
            }
        }
    }
    
    if(lastInput) {
        size_t bit = sizeof(uint32_t) * 8;
        size_t num = (repeat_list_.size() + bit - 1) / bit;
        done_.resizeExact(num);
        done_.fillZero();
        
#ifndef NDEBUG
        string prev_seq = "";
        for(size_t i = 0; i < repeat_index_.size(); i++) {
            TIndexOffU saBegin = repeat_index_[i];
            TIndexOffU saEnd = (i + 1 < repeat_index_.size() ? repeat_index_[i+1] : repeat_list_.size());
            string seq = getString(s, repeat_list_[saBegin], seed_len());
            for(size_t j = saBegin + 1; j < saEnd; j++) {
                string tmp_seq = getString(s, repeat_list_[j], seed_len());
                assert_eq(seq, tmp_seq);
            }
            if(prev_seq != "" ) {
                assert_lt(prev_seq, seq);
            }
            prev_seq = seq;
        }
#endif
    }
}

template<typename TStr>
Range RB_SubSA::find(const TStr& s,
                     const string& seq) const
{
    TIndexOffU i = find_repeat_idx(s, seq);
    if(i >= repeat_index_.size())
        return Range(0, 0);
    
    Range range(repeat_index_[i],
                i + 1 < repeat_index_.size() ? repeat_index_[i+1] : repeat_list_.size());
    return range;
}

template<typename TStr>
TIndexOffU RB_SubSA::find_repeat_idx(const TStr& s,
                                     const string& seq) const
{
    assert_eq(seq.length(), seed_len_);
    
    string temp;
    size_t l = 0, r = repeat_index_.size();
    while(l < r) {
        size_t m = (r + l) >> 1;
        TIndexOffU saElt = repeat_list_[repeat_index_[m]];
        getString(s, saElt, seed_len_, temp);
        if(seq == temp) {
            return m;
        } else if(seq < temp) {
            r = m;
        } else {
            assert(seq > temp);
            l = m + 1;
        }
    }
    return repeat_index_.size();
}

void RB_SubSA::setDone(TIndexOffU off, TIndexOffU len)
{
    assert_leq(off + len, repeat_list_.size());
    const TIndexOffU bit = sizeof(uint32_t) * 8;
    for(TIndexOffU i = off; i < off + len; i++) {
        TIndexOffU quotient = i / bit;
        TIndexOffU remainder = i % bit;
        assert_lt(quotient, done_.size());
        uint32_t num = done_[quotient];
        num = num | (1 << remainder);
        done_[quotient] = num;
    }
}

bool RB_SubSA::isDone(TIndexOffU off, TIndexOffU len) const
{
    assert_leq(off + len, repeat_list_.size());
    const TIndexOffU bit = sizeof(uint32_t) * 8;
    for(TIndexOffU i = off; i < off + len; i++) {
        TIndexOffU quotient = i / bit;
        TIndexOffU remainder = i % bit;
        assert_lt(quotient, done_.size());
        uint32_t num = done_[quotient];
        num = num & (1 << remainder);
        if(num == 0)
            return false;
    }
    
    return true;
}

template<typename TStr>
void RB_SubSA::buildRepeatBase(const TStr& s,
                               CoordHelper& coordHelper,
                               const size_t max_len,
                               EList<RB_RepeatBase>& repeatBases)
{
    if(repeat_index_.empty())
        return;
    
    done_.fillZero();
    
    EList<uint8_t> senseDominant;
    senseDominant.resizeExact(repeat_index_.size());
    senseDominant.fillZero();
    
    EList<size_t> repeatStack;
    EList<pair<TIndexOffU, TIndexOffU> > size_table;
    size_table.reserveExact(repeat_index_.size() / 2 + 1);
    for(size_t i = 0; i < repeat_index_.size(); i++) {
        TIndexOffU begin = repeat_index_[i];
        TIndexOffU end = (i + 1 < repeat_index_.size() ? repeat_index_[i+1] : repeat_list_.size());
        assert_lt(begin, end);
        EList<TIndexOffU> positions; positions.reserveExact(end - begin);
        for(size_t j = begin; j < end; j++) positions.push_back(repeat_list_[j]);
        
        if(!isSenseDominant(coordHelper, positions, seed_len_))
            continue;
        
        senseDominant[i] = 1;
        size_table.expand();
        size_table.back().first = end - begin;
        size_table.back().second = i;
    }
    size_table.sort();

    size_t bundle_count = 0;
    string tmp_str;
    EList<Range> tmp_ranges; tmp_ranges.resizeExact(4);
    EList<pair<size_t, size_t> > tmp_sort_ranges; tmp_sort_ranges.resizeExact(4);
    for(int64_t i = (int64_t)size_table.size() - 1; i >= 0; i--) {
        TIndexOffU idx = size_table[i].second;
        ASSERT_ONLY(TIndexOffU num = size_table[i].first);
        assert_lt(idx, repeat_index_.size());
#ifndef NDEBUG
        if(idx + 1 < repeat_index_.size()) {
            assert_eq(repeat_index_[idx] + num, repeat_index_[idx+1]);
        } else {
            assert_eq(repeat_index_[idx] + num, repeat_list_.size());
        }
#endif
        TIndexOffU saBegin = repeat_index_[idx];
        if(isDone(saBegin)) {
            assert(isDone(saBegin, num));
            continue;
        }

        ASSERT_ONLY(size_t rb_done = 0);
        assert_lt(saBegin, repeat_list_.size());
        repeatStack.push_back(idx);
        while(!repeatStack.empty()) {
            TIndexOffU idx = repeatStack.back();
            assert(senseDominant[idx]);
            repeatStack.pop_back();
            assert_lt(idx, repeat_index_.size());
            TIndexOffU saBegin = repeat_index_[idx];
            TIndexOffU saEnd = (idx + 1 < repeat_index_.size() ? repeat_index_[idx + 1] : repeat_list_.size());
            if(isDone(saBegin)) {
                assert(isDone(saBegin, saEnd - saBegin));
                continue;
            }
            
            TIndexOffU saElt = repeat_list_[saBegin];
            size_t ri = repeatBases.size();
            repeatBases.expand();
            repeatBases.back().seq = getString(s, saElt, seed_len_);
            repeatBases.back().nodes.clear();
            repeatBases.back().nodes.push_back(idx); 
            ASSERT_ONLY(rb_done++);
            setDone(saBegin, saEnd - saBegin);
            bool left = true;
            while(repeatBases[ri].seq.length() <= max_len) {
                if(left) {
                    tmp_str = "N";
                    tmp_str += repeatBases[ri].seq.substr(0, seed_len_ - 1);
                } else {
                    tmp_str = repeatBases[ri].seq.substr(repeatBases[ri].seq.length() - seed_len_ + 1, seed_len_ - 1);
                    tmp_str.push_back('N');
                }
                assert_eq(tmp_str.length(), seed_len_);
                
                for(size_t c = 0; c < 4; c++) {
                    if(left) tmp_str[0] = "ACGT"[c];
                    else     tmp_str.back() = "ACGT"[c];
                    assert_eq(tmp_str.length(), seed_len_);
                    TIndexOffU idx = find_repeat_idx(s, tmp_str);
                    size_t num = 0;
                    if(idx < repeat_index_.size()) {
                        if(idx + 1 < repeat_index_.size()) {
                            num = repeat_index_[idx+1] - repeat_index_[idx];
                        } else {
                            num = repeat_list_.size() - repeat_index_[idx];
                        }
                    }
                    tmp_ranges[c].first = idx;
                    tmp_ranges[c].second = num;
                    assert(num == 0 || num >= seed_count_);
                    tmp_sort_ranges[c].first = num;
                    tmp_sort_ranges[c].second = c;
                    if(idx == repeat_index_.size() ||
                       isDone(repeat_index_[idx]) ||
                       !senseDominant[idx]) {
#ifndef NDEBUG
                        if(idx < repeat_index_.size()) {
                            assert(isDone(repeat_index_[idx], num) ||
                                   !senseDominant[idx]);
                        }
#endif
                        tmp_sort_ranges[c].first = 0;
                        tmp_ranges[c].second = 0;
                    }
                }
                tmp_sort_ranges.sort();
                if(tmp_sort_ranges[3].first < seed_count_) {
                    if(left) {
                        left = false;
                        continue;
                    } else {
                        break;
                    }
                }
                
                for(size_t cc = 0; cc < 3; cc++) {
                    assert_leq(tmp_sort_ranges[cc].first, tmp_sort_ranges[cc+1].first);
                    if(tmp_sort_ranges[cc].first < seed_count_)
                        continue;
                    
                    size_t c = tmp_sort_ranges[cc].second;
                    repeatStack.push_back(tmp_ranges[c].first);
                }

                size_t c = tmp_sort_ranges[3].second;
                if(repeatBases[ri].seq.length() >= max_len) {
                    assert_eq(repeatBases[ri].seq.length(), max_len);
                    TIndexOffU idx = tmp_ranges[c].first;
                    assert(!isDone(repeat_index_[idx]));
                    repeatStack.push_back(idx);
                    if(left) {
                        left = false;
                        continue;
                    } else {
                        break;
                    }
                } else {
                    TIndexOffU idx = tmp_ranges[c].first;
                    TIndexOffU num = tmp_ranges[c].second;
                    setDone(repeat_index_[idx], num);
                    if(left) {
                        repeatBases[ri].seq.insert(0, 1, "ACGT"[c]);
                        repeatBases[ri].nodes.insert(idx, 0);
                    } else {
                        repeatBases[ri].seq.push_back("ACGT"[c]);
                        repeatBases[ri].nodes.push_back(idx);
                    }
                }
            }
        }

        assert_gt(rb_done, 0);
        bundle_count++;
    }

    cerr << "Bundle count: " << bundle_count << endl;

#ifndef NDEBUG
    {
        set<TIndexOffU> idx_set;
        for(size_t i = 0; i < repeatBases.size(); i++) {
            const RB_RepeatBase& repeatBase = repeatBases[i];
            for(size_t j = 0; j < repeatBase.nodes.size(); j++) {
                assert(idx_set.find(repeatBase.nodes[j]) == idx_set.end());
                idx_set.insert(repeatBase.nodes[j]);
            }
        }
    }
#endif
    
#if 0
    {
        EList<pair<size_t, size_t> > kmer_table;
        string query;
        size_t total = 0, match = 0;
        size_t interval = 1;
        if(repeat_index_.size() >= 10000) {
            interval = repeat_index_.size() / 10000;
        }
        EList<TIndexOffU> positions;
        for(size_t i = 0; i < repeat_index_.size(); i += interval) {
            TIndexOffU saElt_idx = repeat_index_[i];
            TIndexOffU saElt_idx_end = (i + 1 < repeat_index_.size() ? repeat_index_[i+1] : repeat_list_.size());
            positions.clear();
            for(size_t j = saElt_idx; j < saElt_idx_end; j++) {
                positions.push_back(repeat_list_[j]);
#ifndef NDEBUG
                if(j > saElt_idx) {
                    TIndexOffU lcp_len = getLCP(s,
                                                coordHelper,
                                                positions[0],
                                                positions.back(),
                                                seed_len_);
                    assert_geq(lcp_len, seed_len_);
                }
                
                TIndexOffU saElt = repeat_list_[j];
                TIndexOffU start = coordHelper.getStart(saElt);
                TIndexOffU start2 = coordHelper.getStart(saElt + seed_len_ - 1);
                assert_eq(start, start2);
#endif
            }
            
            TIndexOffU saElt = repeat_list_[saElt_idx];
            size_t true_count = saElt_idx_end - saElt_idx;
            getString(s, saElt, seed_len_, query);
            
            total++;
            
            size_t count = 0;
            for(size_t r = 0; r < repeatBases.size(); r++) {
                const RB_RepeatBase& repeat = repeatBases[r];
                int pos = repeat.seq.find(query);
                if(pos != string::npos) {
                    for(size_t j = 0; j < repeat.nodes.size(); j++) {
                        TIndexOffU _node = repeat.nodes[j];
                        TIndexOffU _saElt_idx = repeat_index_[_node];
                        TIndexOffU _saElt_idx_end = (_node + 1 < repeat_index_.size() ? repeat_index_[_node+1] : repeat_list_.size());
                        TIndexOffU _saElt = repeat_list_[_saElt_idx];
                        string seq = getString(s, _saElt, seed_len());
                        if(query == seq) {
                            count += (_saElt_idx_end - _saElt_idx);
                        }
                    }
                }
            }
            
            size_t rc_count = 0;
            string rc_query = reverse_complement(query);
            for(size_t r = 0; r < repeatBases.size(); r++) {
                const RB_RepeatBase& repeat = repeatBases[r];
                int pos = repeat.seq.find(rc_query);
                if(pos != string::npos) {
                    for(size_t j = 0; j < repeat.nodes.size(); j++) {
                        TIndexOffU _node = repeat.nodes[j];
                        TIndexOffU _saElt_idx = repeat_index_[_node];
                        TIndexOffU _saElt_idx_end = (_node + 1 < repeat_index_.size() ? repeat_index_[_node+1] : repeat_list_.size());
                        TIndexOffU _saElt = repeat_list_[_saElt_idx];
                        string seq = getString(s, _saElt, seed_len());
                        if(rc_query == seq) {
                            rc_count += (_saElt_idx_end - _saElt_idx);
                        }
                    }
                }
            }
            
            if(count == true_count || rc_count == true_count) {
                match++;
            } else if(total - match <= 10) {
                cerr << "query: " << query << endl;
                cerr << "true count: " << true_count << endl;
                cerr << "found count: " << count << endl;
                cerr << "rc found count: " << rc_count << endl;
                cerr << endl;
            }
        }
        
        cerr << "RB_SubSA: sanity check: " << match << " passed (out of " << total << ")" << endl << endl;
    }
#endif
}


#define write_fp(x) fp.write((const char *)&(x), sizeof((x)))

void RB_SubSA::writeFile(ofstream& fp)
{
    write_fp(sa_size_);
    write_fp(seed_len_);
    write_fp(seed_count_);

    size_t sz = temp_suffixes_.size();
    write_fp(sz);
    for(size_t i = 0; i < sz; i++) {
        write_fp(temp_suffixes_[i]);
    }

    sz = repeat_index_.size();
    write_fp(sz);
    for(size_t i = 0; i < sz; i++) {
        write_fp(repeat_index_[i]);
    }

    sz = done_.size();
    write_fp(sz);
    for(size_t i = 0; i < sz; i++) {
        write_fp(done_[i]);
    }
}

#define read_fp(x) fp.read((char *)&(x), sizeof((x)))

void RB_SubSA::readFile(ifstream& fp)
{
    TIndexOffU val;

    read_fp(val);
    rt_assert_eq(val, sa_size_);

    read_fp(val);
    rt_assert_eq(val, seed_len_);

    read_fp(val);
    rt_assert_eq(val, seed_count_);

    size_t sz;
    read_fp(sz);
    temp_suffixes_.resizeExact(sz);
    for(size_t i = 0; i < sz; i++) {
        read_fp(temp_suffixes_[i]);
    }

    size_t val_sz;

    // repeat_index_
    read_fp(val_sz);
    repeat_index_.resizeExact(val_sz);
    for(size_t i = 0; i < val_sz; i++) {
        read_fp(repeat_index_[i]);
    }

    // done
    read_fp(val_sz);
    done_.resizeExact(val_sz);
    for(size_t i = 0; i < val_sz; i++) {
        read_fp(done_[i]);
    }

}


/****************************/
template class RepeatBuilder<SString<char> >;
template void dump_tstr(const SString<char>& );
template bool compareRepeatCoordByJoinedOff(const RepeatCoord<TIndexOffU>& , const RepeatCoord<TIndexOffU>&);
