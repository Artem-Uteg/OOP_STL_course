#include<iostream>
#include<vector>
#include<string>

const int64_t base = 1e9;
const int base_len = 9;

class BigInteger {

private:

    bool is_positive = true;
    std::vector<int64_t> digits;

    void RemoveZeroes();
    void Shift();

    BigInteger Abs() const;

public:

    BigInteger(int64_t big = 0);

    explicit BigInteger(std::string& s);

    bool positive() const;
    size_t digit_size() const;

    int64_t& operator[](int64_t index);
    const int64_t& operator[](int64_t index) const;

    void swap(BigInteger& big);

    static BigInteger GCD(BigInteger left, BigInteger right);

    BigInteger& operator++();
    BigInteger& operator--();
    BigInteger operator++(int);
    BigInteger operator--(int);

    BigInteger& operator%=(const BigInteger& big);
    BigInteger& operator+=(BigInteger big);
    BigInteger& operator*=(const BigInteger& big);
    BigInteger& operator-=(const BigInteger& big);
    BigInteger& operator/=(const BigInteger& big);

    explicit operator bool() const;
    std::string toString() const;
    BigInteger operator-() const;

};

BigInteger::BigInteger(int64_t big) {

    is_positive = (big >= 0);

    if (big == 0) {
        digits.push_back(0);
        return;
    }

    if (!is_positive) big *= -1;

    while (big != 0) {
        digits.push_back(big % base);
        big /= base;
    }
}

BigInteger::BigInteger(std::string& s) {

    if (s.size() == 0) {
        *this = BigInteger();
        return;
    }

    is_positive = (s[0] != '-');
    if (s[0] == '-') s.erase(0, 1);

    int64_t digit;

    while (s.size() > base_len) {
        digit = stoll(s.substr(s.size() - base_len));
        digits.push_back(digit);
        s.erase(s.size() - base_len, s.size() - 1);
    }

    digit = stoll(s);
    digits.push_back(digit);

    while (digits.back() == 0 && digits.size() > 1) digits.pop_back();
    if (digits[0] == 0 && digits.size() == 1) is_positive = true;
}

void BigInteger::swap(BigInteger& big) {
    std::swap(is_positive, big.is_positive);
    digits.swap( big.digits);
}

int64_t& BigInteger::operator[](int64_t index) {
    return digits[index];
}

const int64_t& BigInteger::operator[](int64_t index) const {
    return digits[index];
}

bool BigInteger::positive() const {
    return is_positive;
}

size_t BigInteger::digit_size() const {
    return digits.size();
}

void BigInteger::RemoveZeroes() {
    while (digits.back() == 0 && digits.size() - 1 > 0) {
        digits.pop_back();
    }

    if (digits[0] == 0 && digits.size() == 1) is_positive = true;
}

bool operator<(const BigInteger& left, const BigInteger& right) {

    if (left.positive() != right.positive()) return (right.positive());

    if (left.digit_size() != right.digit_size()) return (left.positive() != (left.digit_size() > right.digit_size()));

    for (int64_t i = static_cast<int64_t>(left.digit_size()) - 1; i >= 0; --i) {
        if (left[i] != right[i]) return (left.positive() != (left[i] > right[i]));
    }

    return false;
}

BigInteger BigInteger::operator-() const {
    BigInteger copy = *this;
    copy.is_positive = copy.digits[0] == 0 && copy.digits.size() == 1 ? true : !is_positive;

    return copy;
}

bool operator==(const BigInteger& left, const BigInteger& right) {
    return !(left < right || right < left);
}

bool operator!=(const BigInteger& left, const BigInteger& right) {
    return !(left == right);
}

bool operator>(const BigInteger& left, const BigInteger& right) {
    return right < left;
}

bool operator<=(const BigInteger& left, const BigInteger& right) {
    return left < right || left == right;
}

bool operator>=(const BigInteger& left, const BigInteger& right) {
    return right < left || left == right;
}

BigInteger& BigInteger::operator+=(BigInteger big) {
    int sgn = 1;

    if (is_positive != big.is_positive) {

        if (!is_positive) {
            if (-big < *this) {
                is_positive = true;
                swap(big);
            }

        } else if (*this < -big) swap(big);

        sgn = -1;
    }

    int64_t carry = 0;
    int64_t add;
    size_t i = 0;

    for (; i < std::min(big.digits.size(), digits.size()); ++i) {

        add = digits[i] + carry + sgn * big.digits[i];
        digits[i] = add % base;
        carry = add / base;

        if (digits[i] >= 0) continue;

        digits[i] += base;
        --carry;

    }

    for (; i < digits.size() && carry != 0; ++i) {
        add = digits[i] + carry;
        digits[i] = add % base;
        carry = add / base;

        if (add >= 0) continue;

        digits[i] += base;
        --carry;
    }

    for (; i < big.digits.size(); ++i) {
        add = big.digits[i] + carry;
        digits.push_back(add % base);
        carry = add / base;
    }

    if (carry > 0) digits.push_back(carry);

    while (digits.back() == 0 && digits.size() > 1) digits.pop_back();
    if (digits[0] == 0 && digits.size() == 1) is_positive = true;

    return *this;
}

BigInteger& BigInteger::operator-=(const BigInteger& big) {
    return *this += -big;
}

BigInteger& BigInteger::operator*=(const BigInteger& big) {
    if (*this == 0 || big == 0) {
        *this = 0;
        return *this;
    }

    bool result_is_positive = (is_positive == big.is_positive);
    std::vector<int64_t> result(digits.size() + big.digits.size() + 1, 0);

    int64_t carry = 0;

    for (size_t i = 0; i < digits.size(); i++) {

        for (size_t j = 0; j < big.digits.size(); ++j) {
            result[i + j] += digits[i] * big.digits[j] + carry;

            carry = result[i + j] / base;
            result[i + j] %= base;
        }

        result[big.digits.size() + i] += carry;
        carry = 0;
    }

    digits = result;
    is_positive = result_is_positive;

    while (digits.size() - 1 > 0 && digits.back() == 0) digits.pop_back();
    if (digits.size() - 1 == 0 && digits[0] == 0) is_positive = true;
    return *this;
}

BigInteger operator+(const BigInteger& left, const BigInteger& right) {
    BigInteger copy = right;
    copy += left;
    return copy;
}

BigInteger operator-(const BigInteger& left, const BigInteger& right) {
    BigInteger copy = left;
    copy -= right;
    return copy;
}

BigInteger operator*(const BigInteger& left, const BigInteger& right) {
    BigInteger copy = left;
    copy *= right;

    return copy;
}

void BigInteger::Shift() {
    if (digits.size() == 0) digits.push_back(0);

    if ((digits[0] == 0) && (digits.size() == 1)) return;

    digits.push_back(0);

    for (size_t i = digits.size() - 1; i > 0; --i) {
        digits[i] = digits[i - 1];
        if (i == 1) digits[0] = 0;
    }
}


BigInteger& BigInteger::operator/=(const BigInteger& big) {

    if (big.Abs() > Abs()) {
        *this = 0;
        return *this;
    }

    BigInteger result, current, temp;
    BigInteger bi = big;
    bi.is_positive = true;

    result.digits.resize(digits.size());

    int64_t x = 0;
    int64_t left;
    int64_t right;

    for (ssize_t i = digits.size() - 1; i >= 0; --i) {

        x = 0;
        current.Shift();

        current.digits[0] = digits[i];
        current.RemoveZeroes();

        int64_t bi_sz = bi.digits.size();
        int64_t cur_sz = current.digits.size();

        int64_t first_digit = (bi_sz < cur_sz) ? current.digits[bi_sz] : 0;
        int64_t second_digit = (bi_sz - 1 < cur_sz) ? current.digits[bi_sz - 1] : 0;

        left = (base * first_digit + second_digit) / (bi.digits.back() + 1);
        right = (base * first_digit + second_digit + 1) / bi.digits.back();

        int64_t middle;
        while (left - right < 1) {

            middle = (left + right) / 2;
            temp = bi * middle;

            if (temp > current)
                right = middle - 1;
            else {
                x = middle;
                left = middle + 1;
            }
        }

        result.digits[i] = x;
        current = current - bi * result.digits[i];
    }

    result.is_positive = (is_positive == big.is_positive);
    result.RemoveZeroes();

    *this = result;
    return *this;
}

BigInteger& BigInteger::operator%=(const BigInteger& big) {
    if (big.Abs() > Abs())
        return *this;

    BigInteger result;
    BigInteger current;
    BigInteger temp;

    BigInteger bi = big;
    bi.is_positive = true;

    result.digits.resize(digits.size());

    int64_t x;
    int64_t left;
    int64_t right;

    for (ssize_t i = digits.size() - 1; i >= 0; --i) {

        current.Shift();
        current.RemoveZeroes();
        current.digits[0] = digits[i];

        int64_t bi_sz = bi.digits.size();
        int64_t cur_sz = current.digits.size();
        x = 0;

        int64_t first_digit = (bi_sz < cur_sz) ? current.digits[bi_sz] : 0;
        int64_t second_digit = (bi_sz - 1 < cur_sz) ? current.digits[bi_sz - 1] : 0;

        left = (base * first_digit + second_digit) / (bi.digits.back() + 1);
        right = (base * first_digit + second_digit + 1) / bi.digits.back();

        int64_t middle;

        while (left - right < 1) {
            middle = (left + right) / 2;
            temp = bi * middle;
            if (temp > current)
                right = middle - 1;
            else {
                x = middle;
                left = middle + 1;
            }
        }
        result.digits[i] = x;
        current = current - bi * result.digits[i];
    }

    current.RemoveZeroes();
    current.is_positive = is_positive;

    *this = current;
    return *this;
}


BigInteger& BigInteger::operator++() {
    *this += 1;
    return *this;
}

BigInteger& BigInteger::operator--() {
    *this -= 1;
    return *this;
}

BigInteger BigInteger::operator++(int) {
    BigInteger copy = *this;
    ++*this;
    return copy;
}

BigInteger BigInteger::operator--(int) {
    BigInteger copy = *this;
    --*this;
    return copy;
}

BigInteger::operator bool() const {
    return *this != 0;
}

std::string BigInteger::toString() const {
    std::string s, add;

    int64_t length = digits.size();

    if (!is_positive) s = "-";

    s += std::to_string(digits.back());

    for (int64_t i = length - 2; i > -1; --i) {
        add = std::to_string(digits[i]);
        s += std::string(base_len - add.size(), '0');
        s += add;
    }

    return s;
}

BigInteger operator/(const BigInteger& left, const BigInteger& right) {
    BigInteger copy = left;
    copy /= right;
    return copy;
}

BigInteger operator%(const BigInteger& left, const BigInteger& right) {
    BigInteger copy = left;
    copy %= right;
    return copy;
}

std::istream &operator>>(std::istream& in, BigInteger& big) {
    std::string s;
    in >> s;
    big = BigInteger(s);
    return in;
}

std::ostream &operator<<(std::ostream& out, const BigInteger& big) {
    return (out << big.toString());
}

BigInteger BigInteger::GCD(BigInteger left, BigInteger right) {
    right.is_positive = left.is_positive = true;

    while (right > 0 && left > 0){
        if (left > right) left %= right;
        else right %= left;

    }

    return left + right;
}

BigInteger BigInteger::Abs() const {
    BigInteger copy = *this;
    copy.is_positive = true;
    return copy;
}

class Rational {

    friend bool operator<(const Rational& left, const Rational& right);

private:

    BigInteger numerator;
    BigInteger denominator;
    void reduce();

public:

    Rational(const BigInteger& n = 0, const BigInteger& d = 1);
    Rational(int64_t n) : Rational(BigInteger(n)) {};

    Rational& operator=(const Rational& rational);
    Rational operator-() const;
    explicit operator double() const;

    Rational& operator+=(const Rational& rational);
    Rational& operator*=(const Rational& rational);
    Rational& operator-=(const Rational& rational);
    Rational& operator/=(const Rational& rational);

    std::string toString() const;
    std::string asDecimal(size_t precision = 0) const;
};

void Rational::reduce() {
    BigInteger x = BigInteger::GCD(numerator, denominator);
    numerator /= x;
    denominator /= x;

    if (denominator < 0) {
        denominator = -denominator;
        numerator = -numerator;
    }
}

Rational::Rational(const BigInteger& n, const BigInteger& d) : numerator(n), denominator(d) {
    reduce();
}

std::string Rational::toString() const {

    std::string s = numerator.toString();
    if (denominator > 1) s += "/" + denominator.toString();

    return s;
}

Rational Rational::operator-() const {
    Rational copy = *this;
    copy.numerator *= -1;
    return copy;
}

Rational& Rational::operator+=(const Rational& rational) {
    numerator = numerator * rational.denominator + rational.numerator * denominator;
    denominator *= rational.denominator;
    reduce();

    return *this;
}

Rational& Rational::operator*=(const Rational& rational) {
    denominator *= rational.denominator;
    numerator *= rational.numerator;
    reduce();
    return *this;
}

Rational& Rational::operator-=(const Rational& rational) {
    numerator = numerator * rational.denominator - rational.numerator * denominator;
    denominator *= rational.denominator;
    reduce();
    return *this;
}

Rational& Rational::operator/=(const Rational& rational) {
    numerator *= rational.denominator;
    denominator *= rational.numerator;
    reduce();
    return *this;
}

Rational operator+(const Rational& left, const Rational& right) {
    Rational copy = left;
    copy += right;
    return copy;
}

Rational operator*(const Rational& left, const Rational& right) {
    Rational copy = left;
    copy *= right;
    return copy;
}

Rational operator-(const Rational& left, const Rational& right) {
    Rational copy = left;
    copy -= right;
    return copy;
}

Rational operator/(const Rational& left, const Rational& right) {
    Rational copy = left;
    copy /= right;
    return copy;
}

std::string Rational::asDecimal(size_t precision) const {

    BigInteger copy_num = numerator;
    BigInteger before_dot = copy_num / denominator;

    std::string ans;

    ans = (before_dot == 0 && copy_num < 0 ? "-" : "");
    ans += before_dot.toString();

    if (copy_num < 0) copy_num *= -1;

    copy_num %= denominator;

    if (precision == 0) return ans;

    ans += ".";
    copy_num *= 10;

    for (size_t i = 0; i < precision; ++i) {
        ans += (copy_num / denominator).toString();
        copy_num %= denominator;
        copy_num *= 10;
    }

    return ans;
}

bool operator<(const Rational& left, const Rational& right) {
    return left.numerator * right. denominator < left.denominator * right.numerator;
}

bool operator>(const Rational& left, const Rational& right) {
    return right < left;
}

bool operator==(const Rational& left, const Rational& right) {
    return !(right < left || left < right);
}

bool operator!=(const Rational& left, const Rational& right) {
    return !(left == right);
}

bool operator<=(const Rational& left, const Rational& right) {
    return left < right || left == right;
}

bool operator>=(const Rational& left, const Rational& right) {
    return right < left || left == right;
}

Rational& Rational::operator=(const Rational& rational) {
    if (this == &rational) return *this;

    denominator = rational.denominator;
    numerator = rational.numerator;
    return *this;
}

Rational::operator double() const {
    return std::stod(asDecimal(32));
}
