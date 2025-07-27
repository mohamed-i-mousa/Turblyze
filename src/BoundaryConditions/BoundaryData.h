#ifndef BOUNDARYDATA_H
#define BOUNDARYDATA_H

#include <string>

#include "Scalar.h"
#include "Vector.h"

enum class BCType {
    FIXED_VALUE,
    FIXED_GRADIENT,
    ZERO_GRADIENT,
    NO_SLIP,
    UNDEFINED
};

enum class BCValueType {
    SCALAR,
    VECTOR,
    UNDEFINED
};

struct BoundaryData {
    
    BCType type = BCType::UNDEFINED;

    // ----- Value storage ----- //
    BCValueType valueType = BCValueType::UNDEFINED;
    Scalar scalarValue = S(0.0);
    Vector vectorValue; 

    // ----- Gradient storage ----- //
    BCValueType gradientType = BCValueType::UNDEFINED;
    Scalar scalarGradient = S(0.0);
    Vector vectorGradient;

    BoundaryData() = default;

    // ----- Setting Methods ----- // 

    void setFixedValue(Scalar s_val) {
        type = BCType::FIXED_VALUE;
        scalarValue = s_val; 
        valueType = BCValueType::SCALAR;

        // Clear other value types
        vectorValue = Vector();
        gradientType = BCValueType::UNDEFINED;
    }

    void setFixedValue(const Vector& v_val) {
        type = BCType::FIXED_VALUE;
        vectorValue = v_val; 
        valueType = BCValueType::VECTOR;

        // Clear other value types
        scalarValue = S(0.0);
        gradientType = BCValueType::UNDEFINED;
    }

    void setFixedGradient(Scalar s_grad) {
        type = BCType::FIXED_GRADIENT;
        scalarGradient = s_grad; 
        gradientType = BCValueType::SCALAR;

        // Clear other value types
        vectorValue = Vector();
        valueType = BCValueType::UNDEFINED;
    }

    void setFixedGradient(const Vector& v_grad) {
        type = BCType::FIXED_GRADIENT;
        vectorGradient = v_grad; 
        gradientType = BCValueType::VECTOR;

        // Clear other value types
        scalarValue = S(0.0);
        valueType = BCValueType::UNDEFINED;
    }

    void setZeroGradient() {
        type = BCType::ZERO_GRADIENT;
        scalarGradient = S(0.0);
        vectorGradient = Vector(S(0.0), S(0.0), S(0.0));

        valueType = BCValueType::UNDEFINED;
        gradientType = BCValueType::UNDEFINED;
    }

    void setNoSlip() {
        type = BCType::NO_SLIP;
        vectorValue = Vector(S(0.0), S(0.0), S(0.0));
        valueType = BCValueType::VECTOR;

        scalarValue = S(0.0);
        gradientType = BCValueType::UNDEFINED;
    }

    // ----- Getting Methods ----- // 

    Scalar getFixedScalarValue() const {
        if (type == BCType::FIXED_VALUE && valueType == BCValueType::SCALAR) {
            return scalarValue;
        }
        throw std::runtime_error("Attempted to get fixed scalar value, but BC is not set to FIXED_VALUE with SCALAR type.");
    }

    const Vector& getFixedVectorValue() const {
        if (type == BCType::FIXED_VALUE && valueType == BCValueType::VECTOR) {
            return vectorValue;
        }
        // No slip is a special case of vector value
        if (type == BCType::NO_SLIP && valueType == BCValueType::VECTOR) {
            return vectorValue; // which should be (0,0,0)
        }
        throw std::runtime_error("Attempted to get fixed vector value, but BC is not set to FIXED_VALUE/NO_SLIP with VECTOR type.");
    }

    Scalar getFixedScalarGradient() const {
        if (type == BCType::FIXED_GRADIENT && gradientType == BCValueType::SCALAR) {
            return scalarGradient;
        }
        throw std::runtime_error("Attempted to get fixed scalar gradient, but BC is not set to FIXED_GRADIENT with SCALAR type.");
    }

    const Vector& getFixedVectorGradient() const {
        if (type == BCType::FIXED_GRADIENT && gradientType == BCValueType::VECTOR) {
            return vectorGradient;
        }
        throw std::runtime_error("Attempted to get fixed vector gradient, but BC is not set to FIXED_GRADIENT with VECTOR type.");
    }
};

#endif