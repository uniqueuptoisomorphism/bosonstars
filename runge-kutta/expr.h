#ifndef EXPR_H_
#define EXPR_H_

#include <iostream>

#include <gmp.h>

class Expr
{
public:
	Expr* left;
	Expr* right;
	size_t refCount;

	Expr(Expr* left, Expr* right) : left(left), right(right), refCount(0)
	{
		if (left)
			left->grab();
		if (right)
			right->grab();
	}

	virtual ~Expr()
	{
	}

	virtual mpf_t* eval(mpf_t vars[]) = 0;

	void grab()
	{
		++refCount;
	}

	virtual void release()
	{
		--refCount;
		if (refCount <= 0)
		{
			if (left)
				left->release();
			if (right)
				right->release();
			delete this;
		}
	}
};

class DataExpr : public Expr
{
public:
	mpf_t temp;

	DataExpr(Expr* left, Expr* right) : Expr(left, right)
	{
		mpf_init(temp);
	}

	virtual ~DataExpr() override
	{
		mpf_clear(temp);
	}
};

typedef void (*mpf_func)(mpf_t, const mpf_t, const mpf_t);
typedef void (*mpf_func_bin)(mpf_t, const mpf_t);

template <typename T>
using mpf_func_cmp = int (*)(const mpf_t, T);

template <typename T>
using mpf_func_set = void (*)(mpf_t, T);

template <mpf_func f>
class ExprFunc : public DataExpr
{
public:
	ExprFunc(Expr* left, Expr* right) : DataExpr(left, right) { }

	virtual mpf_t* eval(mpf_t vars[]) override
	{
		mpf_t* f1 = left->eval(vars);
		mpf_t* f2 = right->eval(vars);
		f(temp, *f1, *f2);
		return &temp;
	}
};

template <mpf_func_bin f>
class ExprFuncBin : public DataExpr
{
public:
	ExprFuncBin(Expr* left) : DataExpr(left, nullptr) { }

	virtual mpf_t* eval(mpf_t vars[]) override
	{
		mpf_t* f1 = left->eval(vars);
		f(temp, *f1);
		return &temp;
	}
};

class ExprDataRoot : public DataExpr
{
public:
	ExprDataRoot() : DataExpr(nullptr, nullptr) { }

	virtual mpf_t* eval(mpf_t vars[]) override
	{
		return &temp;
	}
};

template <typename T>
class ExprDataRootT : public ExprDataRoot { };

template <>
class ExprDataRootT<long int> : public ExprDataRoot
{
public:
	ExprDataRootT(const long int val) { mpf_set_si(temp, val); }
};

template <>
class ExprDataRootT<unsigned long int> : public ExprDataRoot
{
public:
	ExprDataRootT(const unsigned long int val) { mpf_set_ui(temp, val); }
};

template <>
class ExprDataRootT<double> : public ExprDataRoot
{
public:
	ExprDataRootT(const double val) { mpf_set_d(temp, val); }
};

template <>
class ExprDataRootT<int> : public ExprDataRoot
{
public:
	ExprDataRootT(const int val) { mpf_set_si(temp, val); }
};

template <>
class ExprDataRootT<mpf_t> : public ExprDataRoot
{
public:
	ExprDataRootT(const mpf_t val) { mpf_set(temp, val); }
};

class ExprVar : public Expr
{
public:
	size_t varId;

	ExprVar(const size_t varId) : Expr(nullptr, nullptr), varId(varId) { }

	virtual mpf_t* eval(mpf_t vars[]) override
	{
		return &vars[varId];
	}
};

struct ExprPtr
{
public:
	Expr* p;

	ExprPtr(Expr* p) : p(p) { }
	ExprPtr(const long unsigned int v) : p(new ExprDataRootT<long unsigned int>(v)) { }
	ExprPtr(const long int v) : p(new ExprDataRootT<long int>(v)) { }
	ExprPtr(const int v) : p(new ExprDataRootT<long int>((long int)v)) { }
};

template <typename T, mpf_func_cmp<T> f>
class ExprVarConditional : public Expr
{
public:
	size_t varId;
	T value;
	Expr* eqThreshold;

	ExprVarConditional(size_t varId, T value, ExprPtr greaterThreshold, ExprPtr eqThreshold, ExprPtr lessThreshold) :
		Expr(greaterThreshold.p, lessThreshold.p), varId(varId), value(value), eqThreshold(eqThreshold.p)
	{
		this->eqThreshold->grab();
	}

	virtual mpf_t* eval(mpf_t vars[]) override
	{
		auto result = f(vars[varId], value);

		if (result > 0)
			return left->eval(vars);
		else if (result == 0)
			return eqThreshold->eval(vars);
		else
			return right->eval(vars);
	}

	virtual void release() override
	{
		eqThreshold->release();
		Expr::release();
	}
};

inline ExprPtr operator+(ExprPtr ex1, ExprPtr ex2)
{
	return { new ExprFunc<mpf_add>(ex1.p, ex2.p) };
}

inline ExprPtr operator*(ExprPtr ex1, ExprPtr ex2)
{
	return { new ExprFunc<mpf_mul>(ex1.p, ex2.p) };
}

inline ExprPtr operator/(ExprPtr ex1, ExprPtr ex2)
{
	return { new ExprFunc<mpf_div>(ex1.p, ex2.p) };
}

inline ExprPtr operator-(ExprPtr ex1, ExprPtr ex2)
{
	return { new ExprFunc<mpf_sub>(ex1.p, ex2.p) };
}

inline ExprPtr operator-(ExprPtr ex)
{
	return { new ExprFuncBin<mpf_neg>(ex.p) };
}

template <typename T>
ExprPtr operator+(ExprPtr ex1, T ex2)
{
	return { new ExprFunc<mpf_add>(ex1.p, new ExprDataRootT<T>(ex2)) };
}

template <typename T>
ExprPtr operator*(ExprPtr ex1, T ex2)
{
	return { new ExprFunc<mpf_mul>(ex1.p, new ExprDataRootT<T>(ex2)) };
}

template <typename T>
ExprPtr operator/(ExprPtr ex1, T ex2)
{
	return { new ExprFunc<mpf_div>(ex1.p, new ExprDataRootT<T>(ex2)) };
}

template <typename T>
ExprPtr operator-(ExprPtr ex1, T ex2)
{
	return { new ExprFunc<mpf_sub>(ex1.p, new ExprDataRootT<T>(ex2)) };
}

template <typename T>
ExprPtr operator+(T ex1, ExprPtr ex2)
{
	return { new ExprFunc<mpf_add>(new ExprDataRootT<T>(ex1), ex2.p) };
}

template <typename T>
ExprPtr operator*(T ex1, ExprPtr ex2)
{
	return { new ExprFunc<mpf_mul>(new ExprDataRootT<T>(ex1), ex2.p) };
}

template <typename T>
ExprPtr operator/(T ex1, ExprPtr ex2)
{
	return { new ExprFunc<mpf_div>(new ExprDataRootT<T>(ex1), ex2.p) };
}

template <typename T>
ExprPtr operator-(T ex1, ExprPtr ex2)
{
	return { new ExprFunc<mpf_sub>(new ExprDataRootT<T>(ex1), ex2.p) };
}

inline ExprPtr sqr(const ExprPtr x) { return x * x; }
template <size_t n>
inline ExprPtr ipow(const ExprPtr x) { return x * ipow<n - 1>(x); }
template <>
inline ExprPtr ipow<1>(const ExprPtr x) { return x; }

static ExprPtr Power(const ExprPtr x, const int n) { return (n == 1) ? x : x * Power(x, n - 1); }

#endif
