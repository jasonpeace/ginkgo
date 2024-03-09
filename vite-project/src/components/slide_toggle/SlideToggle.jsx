import {useState, useEffect} from 'react';
import useDebounce from '../use_debounce/UseDebounce'

const SlideToggle = ({toggleOn=false, handleUpdate, isRound}) => {


const outer_style={
    height:"50px",
    width:"100px",
    padding: "5px",
    backgroundColor:"blue",
    position: "relative",
};

const inner_style={
    top:"5px",
    backgroundColor:"red",
    height:"40px",
    width:"40px",
    position: "relative",
    transition: "0.5s",
};

const styleOff={
    ...inner_style,
    left:"5px",
};

const styleOn={
    ...inner_style,
    left:"55px",
};


const getScrubberStyles = (isToggleOn) =>{
    if (isToggleOn){
        var derivedStyle =  styleOn;
    } else{
        var derivedStyle =  styleOff;
    }
    if (isRound){
        derivedStyle["borderRadius"] = "50px"
    }
    return derivedStyle;
}

const getTrackStyles = (isToggleOn) =>{
    let derivedStyle = {...outer_style};
    if (isRound){
        derivedStyle["borderRadius"] = "50px"
    }
    return derivedStyle;
}


const [scrubberStyle, setScrubberStyle] = useState(getScrubberStyles(toggleOn))
const debouncedScrubberStyle = useDebounce(scrubberStyle, 400)
const [trackStyle, setTrackStyle] = useState(getTrackStyles())


const handleClick = () => {
    const new_val = !toggleOn;
    handleUpdate(new_val)

    setScrubberStyle(getScrubberStyles(new_val))
}

return(
    <div onClick={handleClick} style={trackStyle}>
        <div style={debouncedScrubberStyle}/>
    </div>
    );

}

export default SlideToggle
