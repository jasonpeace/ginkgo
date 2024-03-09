import { useEffect, useState } from 'react';

const useDebounce = (value, delay=500) => {
    const [debouncedValue, setDebouncedValue] = useState(value)

    useEffect(()=>{
        const t = setTimeout(()=>{
            setDebouncedValue(value)
        }, delay);

        return() => clearTimeout(t)

    },[value, delay]);

    return debouncedValue

}

export default useDebounce
