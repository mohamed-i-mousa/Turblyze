#include "BoundaryData.h"
#include <stdexcept>

void BoundaryData::setFixedValue(Scalar s_val) 
{
    type = BCType::FIXED_VALUE;
    scalarValue = s_val; 
    valueType = BCValueType::SCALAR;

    vectorValue = Vector();
    gradientType = BCValueType::UNDEFINED;
}

void BoundaryData::setFixedValue(const Vector& v_val) 
{
    type = BCType::FIXED_VALUE;
    vectorValue = v_val; 
    valueType = BCValueType::VECTOR;

    scalarValue = S(0.0);
    gradientType = BCValueType::UNDEFINED;
}

void BoundaryData::setFixedGradient(Scalar s_grad) 
{
    type = BCType::FIXED_GRADIENT;
    scalarGradient = s_grad; 
    gradientType = BCValueType::SCALAR;

    vectorValue = Vector();
    valueType = BCValueType::UNDEFINED;
}

void BoundaryData::setFixedGradient(const Vector& v_grad) 
{
    type = BCType::FIXED_GRADIENT;
    vectorGradient = v_grad; 
    gradientType = BCValueType::VECTOR;

    scalarValue = S(0.0);
    valueType = BCValueType::UNDEFINED;
}

void BoundaryData::setZeroGradient() 
{
    type = BCType::ZERO_GRADIENT;
    scalarGradient = S(0.0);
    vectorGradient = Vector(S(0.0), S(0.0), S(0.0));

    valueType = BCValueType::UNDEFINED;
    gradientType = BCValueType::UNDEFINED;
}


void BoundaryData::setNoSlip() 
{
    type = BCType::NO_SLIP;
    vectorValue = Vector(S(0.0), S(0.0), S(0.0));
    valueType = BCValueType::VECTOR;

    scalarValue = S(0.0);
    gradientType = BCValueType::UNDEFINED;
}

Scalar BoundaryData::getFixedScalarValue() const 
{
    if (type == BCType::FIXED_VALUE && valueType == BCValueType::SCALAR) 
    {
        return scalarValue;
    }
    
    throw std::runtime_error
    (
        "Attempted to get fixed scalar value, but BC is not set to "
        "FIXED_VALUE with SCALAR type."
    );
}

const Vector& BoundaryData::getFixedVectorValue() const 
{
    if (type == BCType::FIXED_VALUE && valueType == BCValueType::VECTOR)
    {
        return vectorValue;
    }

    if (type == BCType::NO_SLIP && valueType == BCValueType::VECTOR) 
    {
        return vectorValue;
    }

    throw std::runtime_error
    (
        "Attempted to get fixed vector value, but BC is not set to "
        "FIXED_VALUE/NO_SLIP with VECTOR type."
    );
}

Scalar BoundaryData::getFixedScalarGradient() const 
{
    if 
    (
        type == BCType::FIXED_GRADIENT 
     && gradientType == BCValueType::SCALAR
    ) 
    {
        return scalarGradient;
    }

    throw std::runtime_error
    (
        "Attempted to get fixed scalar gradient, but BC is not set to "
        "FIXED_GRADIENT with SCALAR type."
    );
}

const Vector& BoundaryData::getFixedVectorGradient() const 
{
    if 
    (
        type == BCType::FIXED_GRADIENT 
     && gradientType == BCValueType::VECTOR
    ) 
    {
        return vectorGradient;
    }

    throw std::runtime_error
    (
        "Attempted to get fixed vector gradient, but BC is not set to "
        "FIXED_GRADIENT with VECTOR type."
    );
}
