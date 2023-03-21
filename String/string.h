#include <iostream>
#include <cstring>

class String {
private:
    char* str;
    size_t len;
    size_t memory;

    void modify(size_t size) {
        char* new_str;
        new_str = new char[size];
        memory = size;
        memcpy(new_str, str, len);
        delete[] str;
        str = new_str;
    }
public:
    String(const char* def) : str(new char[strlen(def) * 2]), len(strlen(def)), memory(strlen(def) * 2) {
        memcpy(str, def, len);
    }
    String(const char c) : str(new char{ c }), len(1), memory(1) {}

    String(size_t n, char c) : str(new char[n * 2]), len(n), memory(n * 2) {
        memset(str, c, n);
    }

    String() : str(new char[1]), len(0), memory(1) {};

    size_t length() const {
        return len;
    }

    String(const String& s) : str(new char[s.length() * 2]), len(s.length()), memory(s.length() * 2) {
        memcpy(str, s.str, s.length());
    }
    ~String() {
        if (len != 0) delete[] str;
    }
    void swap(String& other) {
        std::swap(len, other.len);
        std::swap(memory, other.memory);
        std::swap(str, other.str);
    }
    String& operator=(const String& other) {
        if (*this == other) return *this;

        delete[] str;
        len = other.len;
        memory = other.memory;
        str = new char[memory];

        copy(other.str, other.str + len, str);
        return *this
    }


    void pop_back() {
        if (len != 0) {
            --len;
        }

    }

    void push_back(const char elem) {
        if (len == memory) {
            modify(memory * 2);
        }
        ++len;
        str[len - 1] = elem;
    }

    char& front() {
        return str[0];
    }

    const char& front() const {
        return str[0];
    }

    char& back() {
        return str[len - 1];
    }

    const char& back() const {
        return str[len - 1];
    }

    String& operator+=(const String& other) {
        if (len + other.length() >= memory) {
            this->modify((len + other.length()) * 2);
        }

        for (size_t i = len; i < len + other.length(); ++i) {
            str[i] = other[i - len];
        }
        len += other.length();

        return *this;
    }

    String& operator+=(char other) {
        push_back(other);
        return *this;
    }

    char& operator[](size_t index) {
        return str[index];
    }

    char operator[](size_t index) const {
        return str[index];
    }

    size_t find(const String& sub) const {
        for (size_t i = 0; i < len - sub.length() + 1; ++i) {
            if (this->substr(i, sub.length()) == sub)
                return i;
        }
        return len;
    }

    size_t rfind(const String& substr) const {
        for (size_t i = len - substr.length(); static_cast<int>(i) >= 0; --i) {
            bool flag = true;
            for (size_t j = i; j < i + substr.length(); ++j) {
                if (str[j] != substr[j - i]) {
                    flag = false;
                    break;
                }
            }
            if (flag)
                return i;
        }
        return len;
    }

    size_t find(const char* sub) {
        for (size_t i = 0; i < len - strlen(sub) + 1; ++i) {
            if (this->substr(i, strlen(sub)) == sub)
                return i;
        }
        return len;
    }

    size_t rfind(const char* substr) {
        for (size_t i = len - 1; i > (strlen(substr)) - 1; --i) {
            bool flag = true;
            for (size_t j = i; j > i - strlen(substr) + 1; --j) {
                if (str[j] != substr[j - i]) {
                    flag = false;
                    break;
                }
            }
            if (flag)
                return i;
        }
        return len;
    }

    String substr(size_t index, size_t count) const {
        String temp;
        for (size_t i = 0; i < count; ++i) {
            temp.push_back(str[index + i]);
        }
        return temp;
    }

    bool empty() const {
        return len == 0;
    }

    void clear() {
        if (len == 0) return;
        delete[] str;
        str = new char[10];
        memory = 10;
        len = 0;
    }

    bool operator==(const String& other) {
        if (len != other.length())
            return false;
        for (size_t i = 0; i < len; ++i) {
            if (str[i] != other[i]) return false;
        }
        return true;
    }

};

String operator+(const String& first, const String& second) {
    String temp = first;
    temp += second;
    return temp;
}

String operator+(const char& first, const String& second) {
    String temp;
    temp.push_back(first);
    return temp + second;
}

String operator+(const String& first, const char& second) {
    String temp = first;
    temp.push_back(second);
    return temp;
}

std::ostream& operator<<(std::ostream& out, const String& s) {
    for (size_t i = 0; i < s.length(); ++i) {
        out << s[i];
    }
    return out;
}

std::istream& operator>>(std::istream& in, String& s) {
    char* test = new char[1];
    s.clear();
    test[0] = in.get();
    while (!(std::isspace(test[0]) || in.eof())) {
        s.push_back(test[0]);
        test[0] = in.get();
    }

    return in;
}
