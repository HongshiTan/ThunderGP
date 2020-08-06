
>Upon receiving edge __*ei*__ ,__*i*__ ≥ 1:

> if coin(1/__*i*__) = “head” then
>>(__*r1*__, __*r2*__, __*t*__, __*c*__) ← (__*ei*__, __*∅*__, __*∅*__, __*0*__)    // __*e1*__ is the new level-1 edge.

>else

>>if __*ei*__ is adjacent to  __*r1*__ then
>>>__*c*__ ← __*c*__ + 1;

>>>if coin(1/__*c*__) = “head” then
>>>>(__*r2*__ , __*t*__) ← (__*ei*__, __*∅*__) // *__ei__* is the new level-2 edge.

>>>else
>>>>if __*ei*__ forms a triangle with __*r1*__ and __*r2*__ then
>>>>>__*t*__ ← {__*r1*__, __*r2*__, __*ei*__ } // __*ei*__ closes the wedge __*r1*__, __*r2*__ 

>>>>end

>>>end

>>end

>end
