#include <iostream>
#include <fstream>
#include <string>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <omp.h>

using namespace std;

// ---------------------------------------------
// Reverse complement (used only at k-mer load)
// ---------------------------------------------
string reverse_complement(const string& seq) {
    string rc = seq;
    reverse(rc.begin(), rc.end());
    for (char &c : rc) {
        switch (c) {
            case 'A': c = 'T'; break;
            case 'T': c = 'A'; break;
            case 'C': c = 'G'; break;
            case 'G': c = 'C'; break;
            default:  break;   // leave N, etc. unchanged
        }
    }
    return rc;
}

// ---------------------------------------------
// Load kmers into hash set (both strands)
// ---------------------------------------------
unordered_set<string> load_kmers(const string& filename) {
    unordered_set<string> kmers;
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening k-mer file: " << filename << endl;
        return kmers;
    }

    string line;
    while (getline(file, line)) {
        if (line.empty()) continue;
        string rc = reverse_complement(line);
        kmers.insert(line);
        kmers.insert(rc);
    }
    return kmers;
}

// ---------------------------------------------
// Count how many k-mers from the set occur in read
// (reuse a single k-mer buffer, no canonicalization)
// ---------------------------------------------
int count_kmers_in_read(const string& read,
                        const unordered_set<string>& kmers,
                        int klen)
{
    if (read.size() < static_cast<size_t>(klen)) return 0;

    int count = 0;
    string kmer(klen, 'N');  // reusable buffer

    const size_t limit = read.size() - klen;
    for (size_t i = 0; i <= limit; ++i) {
        // copy into reusable kmer buffer
        std::copy(read.begin() + i,
                  read.begin() + i + klen,
                  kmer.begin());

        if (kmers.find(kmer) != kmers.end())
            ++count;
    }
    return count;
}

// ---------------------------------------------
// MAIN
// ---------------------------------------------
int main(int argc, char* argv[]) {

    if (argc != 5) {
        cerr << "Usage: " << argv[0]
             << " <kmer_file> <fastq_file> <output_file> <threads>\n";
        return 1;
    }

    const string kmer_file   = argv[1];
    const string fastq_file  = argv[2];
    const string output_file = argv[3];
    const int threads        = std::atoi(argv[4]);
    const size_t chunk_size  = 1000000;  // number of reads per chunk

    if (threads > 0)
        omp_set_num_threads(threads);

    // Load k-mers
    unordered_set<string> kmers = load_kmers(kmer_file);
    if (kmers.empty()) {
        cerr << "ERROR: No k-mers loaded.\n";
        return 1;
    }

    int klen = kmers.begin()->size();
    cerr << "Loaded " << kmers.size()
         << " k-mers of length " << klen << "\n";

    // Open I/O
    ifstream fastq(fastq_file);
    if (!fastq.is_open()) {
        cerr << "Error: cannot open FASTQ: " << fastq_file << "\n";
        return 1;
    }

    ofstream out(output_file);
    if (!out.is_open()) {
        cerr << "Error: cannot open output file: " << output_file << "\n";
        return 1;
    }

    string header, seq, plus, qual;
    size_t processed = 0;

    vector<string> headers;
    vector<string> seqs;
    headers.reserve(chunk_size);
    seqs.reserve(chunk_size);

    // Read FASTQ in chunks
    while (getline(fastq, header) &&
           getline(fastq, seq)    &&
           getline(fastq, plus)   &&
           getline(fastq, qual))
    {
        headers.push_back(header);
        seqs.push_back(seq);

        // Process when chunk is full
        if (headers.size() >= chunk_size) {

            vector<string> results(headers.size());

            #pragma omp parallel for
            for (size_t i = 0; i < headers.size(); ++i) {
                // remove '@'
                string read_id = (headers[i].empty() || headers[i][0] != '@')
                                 ? headers[i]
                                 : headers[i].substr(1);

                int count = count_kmers_in_read(seqs[i], kmers, klen);
                results[i] = read_id + "\t" + to_string(count);
            }

            for (const auto& s : results)
                out << s << "\n";

            processed += headers.size();
            cerr << "Processed " << processed << " reads...\n";

            headers.clear();
            seqs.clear();
        }
    }

    // Process remaining reads
    if (!headers.empty()) {
        vector<string> results(headers.size());

        #pragma omp parallel for
        for (size_t i = 0; i < headers.size(); ++i) {
            string read_id = (headers[i].empty() || headers[i][0] != '@')
                             ? headers[i]
                             : headers[i].substr(1);

            int count = count_kmers_in_read(seqs[i], kmers, klen);
            results[i] = read_id + "\t" + to_string(count);
        }

        for (const auto& s : results)
            out << s << "\n";

        processed += headers.size();
    }

    cerr << "Finished. Total reads processed: " << processed << "\n";

    fastq.close();
    out.close();
    return 0;
}