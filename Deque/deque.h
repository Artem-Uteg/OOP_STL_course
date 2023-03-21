#include <vector>
#include <stdexcept>

template <typename T>

class Deque {
public:
    Deque();
    Deque(const Deque<T>& other);
    Deque(const int& count);
    Deque(const int& count, const T& elem);

    ~Deque();

    Deque& operator=(const Deque& other);

    size_t size() const { return bucket_size*(back_i - front_i) + back_j - front_j + 1; };

    T& operator[](size_t index);
    const T& operator[](size_t index) const;


    void push_back(const T& elem);
    void pop_back();
    void push_front(const T& elem);

    void pop_front();
    T& at(size_t index);
    const T& at(size_t index) const;

    struct iterator {
    public:

        iterator(int i, int j, Deque<T>* deque) : deque(deque), i(i), j(j)  {} ;
        iterator& operator++() { return *this += 1; };
        iterator& operator--() { return *this -= 1; };
        iterator operator++(int) { return (++*this) - 1; };
        iterator operator--(int) { return (--*this) + 1; };

        iterator& operator+=(const int& n);
        iterator& operator-=(const int& n);
        iterator operator+(const int& n) const;
        iterator operator-(const int& n) const;

        bool operator==(const iterator& it) const { return ((i == it.i) && (j == it.j)); };
        bool operator!=(const iterator& it) const { return !(*this == it); };
        bool operator<(const iterator& it) const { return ((i < it.i) || ((i == it.i) && (j < it.j))); } ;
        bool operator>(const iterator& it) const { return (it < *this); };
        bool operator>=(const iterator& it) const { return ((*this == it) || (*this > it)); };
        bool operator<=(const iterator& it) const { return ((*this == it) || (*this < it)); };

        int operator-(const iterator& it) const { return j + i * bucket_size - (it.i * bucket_size + it.j); };

        T& operator*() { return deque->data[i][j]; };
        T* operator->() { return deque->data[i] + j; };

    private:
        Deque<T>* deque;
        int i = 0;
        int j = 0;
    };

    struct const_iterator {
    public:

        const_iterator(int i, int j, const Deque<T>* deque) : deque(deque), i(i), j(j) {};
        const_iterator(const iterator& it): i(it.i), j(it.j), deque(it.deque) {};

        const_iterator& operator++() { return *this += 1; };
        const_iterator& operator--() { return *this -= 1; };
        const_iterator operator++(int) {return (++*this) - 1; };
        const_iterator operator--(int) {return (--*this) + 1; };
        const_iterator& operator+=(const int& n);
        const_iterator& operator-=(const int& n);
        const_iterator operator+(const int& n) const;
        const_iterator operator-(const int& n) const;

        bool operator==(const const_iterator& it) const { return ((i == it.i) && (j == it.j)); };
        bool operator!=(const const_iterator& it) const { return !(*this == it); } ;
        bool operator<(const const_iterator& it) const { return ((i < it.i) ||((i == it.i) && (j < it.j))); };
        bool operator>(const const_iterator& it) const { return (it < *this);} ;
        bool operator>=(const const_iterator& it) const { return ((*this == it) || (*this > it)); };
        bool operator<=(const const_iterator& it) const { return ((*this == it) || (*this < it)); } ;

        int operator-(const const_iterator& it) const { return j + i * bucket_size - (it.i * bucket_size + it.j); };

        const T& operator*() const {return deque->data[i][j]; } ;
        const T* operator->() const { return deque -> data[i] + j; } ;

    private:
        const Deque<T>* deque;
        int i = 0;
        int j = 0;
    };

    using reverse_iterator = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;
    iterator begin() { return iterator(front_i, front_j, this); };
    const_iterator begin() const { return const_iterator(front_i, front_j, this); };
    iterator end() { return iterator(back_i, back_j, this) + 1; };
    const_iterator end() const { return const_iterator(back_i, back_j, this) + 1; } ;
    const_iterator cbegin() const { return const_iterator(front_i, front_j, this); };
    const_iterator cend() const { return const_iterator(back_i, back_j, this) + 1; };

    void insert(typename Deque<T>::iterator it, const T& elem);
    void erase(typename Deque<T>::iterator it);

    reverse_iterator rbegin() { return reverse_iterator(front_i, front_j, this); };
    reverse_iterator rend() { return reverse_iterator(back_i, back_j, this); };

    const_reverse_iterator crbegin() const{ return const_reverse_iterator(front_i, front_j, this); };
    const_reverse_iterator crend() const { return const_reverse_iterator(back_i, back_j, this); };

private:
    std::vector<T*> data;

    size_t front_i = 1, front_j = 0;
    size_t back_i = 0, back_j = bucket_size - 1;

    static const size_t bucket_size = 32;
};

template <typename T>
Deque<T>::Deque() {
    for (size_t i = 0; i < 2; ++i)
        data.push_back(nullptr);
}

template <typename T>
Deque<T>::Deque(const int& count, const T& elem) {

    for (size_t i = 0; i < 2; ++i)
        data.push_back(nullptr);

    for (int i = 0; i < count; ++i)
        push_back(elem);
}

template <typename T>
Deque<T>::Deque(const int& count) {
    data.push_back(nullptr);
    data.push_back(nullptr);

    T temp;

    for (int i = 0; i < count; ++i) {
        push_back(temp);
    }
}

template <typename T>
Deque<T>::Deque(const Deque<T>& other) {

    size_t i = 0, j = 0;

    back_i = other.back_i;
    back_j = other.back_j;

    front_i = other.front_i;
    front_j = other.front_j;

    T* row = reinterpret_cast<T*>(new int8_t[bucket_size * sizeof(T)]);

    try {
        for (; i < other.data.size(); ++i) {
            if (other.data[i]) {
                for (; j < bucket_size; ++j)
                    if (((i != front_i) || (j + 1 > front_j)) && ((i != back_i) || (j < back_j + 1)))
                        new(j + row) T(other.data[i][j]);
            }
            data.push_back((other.data[i] ? row: nullptr));
        }
    } catch (...) {
        for (size_t di = 0; di < i + 1; ++di) {
            if (!data[di]) continue;
            size_t size = (di < i) ? bucket_size : j;
            for (size_t dj = 0; dj < size; ++dj)
                (data[di] + dj)->~T();

            delete[] reinterpret_cast<int8_t*>(data[i]);
        }
    }
}

template <typename T>
Deque<T>& Deque<T>::operator=(const Deque& other) {

    if (this == &other) return *this;

    for (size_t i = 0; i < data.size(); ++i) {
        if (data[i] == nullptr) continue;

        for (size_t j = 0; j < bucket_size; ++j)
            if (((i != front_i) || (j + 1 > front_j)) && ((i != back_i) || (j < back_j + 1)))
                data[i][j].~T();
        delete[] reinterpret_cast<int8_t*>(data[i]);
    }

    data.clear();

    size_t i = 0, j = 0;

    front_i = other.front_i;
    back_i = other.back_i;

    front_j = other.front_j;
    back_j = other.back_j;

    T* row = reinterpret_cast<T*>(new int8_t[bucket_size * sizeof(T)]);

    try {
        for (; i < other.data.size(); ++i) {
            if (other.data[i] != nullptr) {
                for (; j < bucket_size; ++j)
                    if (((i != front_i) || (j + 1 > front_j)) && ((i != back_i) || (j < back_j + 1)))
                        new(j + row) T(other.data[i][j]);
            }
            data.push_back((other.data[i] != nullptr ? row: nullptr));
        }

    } catch (...) {

        for (size_t di = 0; di < i + 1; ++di) {
            if (data[di] == nullptr) continue;
            size_t size = (di < i) ? bucket_size: j;
            for (size_t dj = 0; dj <  size; ++dj)
                    (data[di] + dj)->~T();

            delete[] reinterpret_cast<int8_t*>(data[i]);
        }
        throw;
    }
    return *this;
}

template <typename T>
Deque<T>::~Deque() {
    for (size_t i = 0; i < data.size(); ++i) {
        if ((i + 1 <= front_i) || (i - 1 >= back_i) || (!data[i])) continue;

        for (size_t j = 0; j < bucket_size; ++j)
            if (((i > front_i) || (j + 1 > front_j)) && ((i < back_i) || (j < back_j + 1)))
                data[i][j].~T();

        delete[] reinterpret_cast<int8_t*>(data[i]);
    }
}

template <typename T>
T& Deque<T>::operator[](size_t index) {
    return data[front_i + ((front_j + index) / bucket_size)][(front_j + index) % bucket_size];
}

template <typename T>
const T& Deque<T>::operator[](size_t index) const {
    return data[front_i + ((front_j + index) / bucket_size)][(front_j + index) % bucket_size];
}

template <typename T>
T& Deque<T>::at(size_t index) {
    if ((index < size()) && (index >= 0)) return (*this)[index];
    throw std::out_of_range("Index out range");
}

template <typename T>
const T& Deque<T>::at(size_t index) const {
    if ((index < size()) && (index >= 0)) return (*this)[index];
    throw std::out_of_range("Index out range");
}

template <typename T>
void Deque<T>::push_back(const T& elem) {
    if (back_j + 1 < bucket_size) {
        new(data[back_i] + (++back_j)) T(elem);
        return;
    }
    if (back_i + 1 < data.size()) {
        back_j = 0;
        data[++back_i] = reinterpret_cast<T*>(new int8_t[bucket_size * sizeof(T)]);
        try {
            new(data[back_i]) T(elem);
        } catch (...) {
            delete[] reinterpret_cast<int8_t*>(data[back_i--]);
            throw;
        }
        return;
    }
    std::vector<T*> new_data(data.size() * 4, nullptr);

    size_t new_size = new_data.size() / 4;
    for (size_t i = 0; i < data.size(); ++i)
        new_data[i + new_size] = data[i];

    data = new_data;

    back_i += new_size;
    front_i += new_size;

    if (back_j + 1 < bucket_size) {
        new(data[back_i] + (++back_j)) T(elem);

    } else if (back_i < data.size() - 1) {
        data[++back_i] = reinterpret_cast<T*>(new int8_t[bucket_size * sizeof(T)]);

        try {
            new(data[back_i]) T(elem);
        } catch (...) {
            delete[] reinterpret_cast<int8_t*>(data[back_i--]);
            throw;
        }
        back_j = 0;
    }
}

template <typename T>
void Deque<T>::push_front(const T& elem) {

    if (front_j > 0) {
        new(data[front_i] + (--front_j)) T(elem);
        return;
    }

    if (front_i > 0) {
        data[--front_i] = reinterpret_cast<T*>(new int8_t[bucket_size * sizeof(T)]);
        try {
            new(bucket_size - 1 + data[front_i]) T(elem);
        } catch (...) {
            delete[] reinterpret_cast<int8_t*>(data[++front_i]);
            throw;
        }
        front_j = bucket_size - 1;
        return;
    }

    std::vector<T*> new_data(data.size() * 4, nullptr);

    size_t new_size = new_data.size() / 4;
    for (size_t i = 0; i < data.size(); ++i)
        new_data[i + new_size] = data[i];

    data = new_data;

    back_i += new_size;
    front_i += new_size;

    if (front_j > 0) {
        new(data[front_i] + (--front_j)) T(elem);
    } else if (front_i > 0) {
        data[--front_i] = reinterpret_cast<T*>(new int8_t[bucket_size * sizeof(T)]);

        try {
            new(bucket_size - 1 + data[front_i]) T(elem);
        } catch (...) {
            delete[] reinterpret_cast<int8_t*>(data[front_i++]);
            throw;
        }
        front_j = bucket_size - 1;
    }
}

template <typename T>
void Deque<T>::pop_back() {

    data[back_i][back_j].~T();

    if (back_j == 0) {
        delete[] reinterpret_cast<int8_t*>(data[back_i--]);
        back_j = bucket_size;
    }
    --back_j;
}

template <typename T>
void Deque<T>::pop_front() {
    data[front_i][front_j].~T();

    if (front_j == bucket_size - 1) {
        front_j = -1;
        delete[] reinterpret_cast<int8_t*>(data[front_i++]);
    }
    ++front_j;
}

template <typename T>
void Deque<T>::erase(Deque<T>::iterator it) {
    if (it == (end() - 1)) {
        pop_back();
    } else {
        (*it).~T();
        new(&(*it)) T(*(it + 1));
        erase(it + 1);
    }
}

template <typename T>
void Deque<T>::insert(Deque<T>::iterator it, const T& elem) {
    if (it == end()) {
        push_back(elem);
    } else {
        insert(it + 1, *it);
        (*it).~T();
        new(&(*it)) T(elem);
    }
}

template <typename T>
typename Deque<T>::iterator Deque<T>::iterator::operator+(const int& n) const {
    Deque<T>::iterator copy = *this;
    return copy += n;
}

template <typename T>
typename Deque<T>::iterator& Deque<T>::iterator::operator+=(const int& n) {
    i += (j + n) / bucket_size;
    j = (j + n) % bucket_size;
    return *this;
}

template <typename T>
typename Deque<T>::iterator Deque<T>::iterator::operator-(const int& n) const {
    Deque<T>::iterator copy = *this;
    return copy -= n;
}

template <typename T>
typename Deque<T>::iterator& Deque<T>::iterator::operator-=(const int& n) {
    i += (j - n) / bucket_size;
    j = (j - n) % bucket_size;
    return *this;
}

template <typename T>
typename Deque<T>::const_iterator Deque<T>::const_iterator::operator+(const int& n) const {
    Deque<T>::const_iterator copy = *this;
    return copy += n;
}


template <typename T>
typename Deque<T>::const_iterator& Deque<T>::const_iterator::operator+=(const int& n) {
    i += (j + n) / bucket_size;
    j = (j + n) % bucket_size;
    return *this;
}

template <typename T>
typename Deque<T>::const_iterator Deque<T>::const_iterator::operator-(const int& n) const {
    Deque<T>::const_iterator copy = *this;
    return copy -= n;
}

template <typename T>
typename Deque<T>::const_iterator& Deque<T>::const_iterator::operator-=(const int& n) {
    i += (j - n) / bucket_size;
    j = (j - n) % bucket_size;
    return *this;
}

