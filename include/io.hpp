#ifndef IO_HPP
#define IO_HPP

#include <fstream>
#include <string>

namespace io {

/*
 * General-purpose load functions
 * Supported types:
 * - all simple C types (NO ARRAYS)
 * - std::vector
 * 
 * Add specializations for other types
 */

template <typename T>
static std::size_t load(std::istream& istrm, T& val)
{
    static_assert(std::is_standard_layout<T>::value);
    if constexpr (std::is_array<T>::value) throw std::domain_error("[Function load] C arrays are not supported");
    istrm.read(reinterpret_cast<char*>(&val), sizeof(T));
    return sizeof(T);
}

template <typename T, typename Allocator>
static std::size_t load(std::istream& istrm, std::vector<T, Allocator>& vec)
{
    std::size_t n;
    load(istrm, n);
    vec.resize(n);
    std::size_t bytes_read = sizeof(n);
    for (auto& v : vec) bytes_read += load(v);
    // istrm.read(reinterpret_cast<char*>(vec.data()), static_cast<std::streamsize>(sizeof(T) * n));
    return bytes_read;
}

//-----------------------------------------------------------------------------------------------------------------------

template <typename T>
static std::size_t store(T const& val, std::ostream& ostrm)
{
    static_assert(std::is_standard_layout<T>::value);
    if constexpr (std::is_array<T>::value) throw std::domain_error("[Function store] C arrays are not supported");
    ostrm.write(reinterpret_cast<char const*>(&val), sizeof(T));
    return sizeof(T);
}

template <typename T, typename Allocator>
static std::size_t store(std::vector<T, Allocator> const& vec, std::ostream& ostrm)
{
    std::size_t n = vec.size();
    std::size_t bytes_written = store(n, ostrm);
    for (auto const& v : vec) bytes_written += store(v);
    // ostrm.write(reinterpret_cast<char const*>(vec.data()), static_cast<std::streamsize>(sizeof(T) * n));
    return bytes_written;
}

//-------------------------------------------------------------------------------------------------------------------------

class loader 
{
    public:
        loader(std::istream& input_stream);

        template <typename T>
        void apply(T& var);

        template <typename T, class Allocator>
        void apply(std::vector<T, Allocator>& vec);

        std::size_t get_byte_size() const {return istrm.tellg();}
        std::size_t get_byte_size_of_simple_types() const {return num_bytes_pods;}
        std::size_t get_byte_size_of_vectors() const {return num_bytes_vecs_of_pods;}

    private:
        std::istream& istrm;
        size_t num_bytes_pods;
        size_t num_bytes_vecs_of_pods;
};

/*
 * The method is intended to be used by classes implementing a visit method.
 * The amount of bytes read is also saved, extending the functionality of the load() functions. 
 */
template <typename T>
void loader::apply(T& var)
{
    if constexpr (std::is_standard_layout<T>::value) {
        auto nr = load(istrm, var);
        num_bytes_pods += nr;
    } else {
        var.visit(*this); // supposing the to-be saved object has a visit() method
    }
}

/**
 * This specialization for std::vector restricts the more general load() to vectors of basic C types.
 * This is to store their size separate from other types.
 */
template <typename T, typename Allocator>
void loader::apply(std::vector<T, Allocator>& vec)
{
    if constexpr (std::is_standard_layout<T>::value) {
        auto nr = load(istrm, vec);
        num_bytes_vecs_of_pods += nr;
    } else {
        std::size_t n;
        apply(n);
        vec.resize(n);
        for (auto& v : vec) apply(v); // Call apply(), not load() since we want to recursively count the number of bytes
        // [[maybe unused]] load(vec);
    }
} 

//-------------------------------------------------------------------------------------------------------------------------------

class saver 
{
    public:
        saver(std::ostream& output_stream);
        
        template <typename T>
        void apply(T& var);

        template <typename T, class Allocator>
        void apply(std::vector<T, Allocator>& vec);

        std::size_t get_byte_size() const {return ostrm.tellp();}

    private:
        std::ostream& ostrm;
};

template <typename T>
void saver::apply(T& var) 
{
    if constexpr (std::is_standard_layout<T>::value) {
        store(var, ostrm);
    } else {
        var.visit(*this);
    }
}

template <typename T, typename Allocator>
void saver::apply(std::vector<T, Allocator>& vec) 
{
    if constexpr (std::is_standard_layout<T>::value) {
        store(vec, ostrm);
    } else {
        size_t n = vec.size();
        apply(n);
        for (auto& v : vec) apply(v);
    }
}

//----------------------------------------------------------------------------------------------------------------------------------

template <typename Visitor, typename T, class StreamType>
static std::size_t visit(T& data_structure, StreamType& strm)
{
    Visitor visitor(strm);
    visitor.apply(data_structure);
    return visitor.get_byte_size();
}

template <typename T>
static std::size_t load(T& data_structure, std::string filename)
{
    std::ifstream strm(filename, std::ios::binary);
    return visit<loader, T>(data_structure, strm);
}

template <typename T>
static std::size_t store(T& data_structure, std::string filename)
{
    std::ofstream strm(filename, std::ios::binary);
    return visit<saver, T>(data_structure, strm);
}

} // namespace io

#endif // IO_HPP