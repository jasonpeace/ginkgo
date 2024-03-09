import { useState, useEffect } from 'react';


export default function Toggle({toggleOn=false, handleUpdate}){
    let on_style = {
        textAlign: "right",
        height:"80px",
        width:"200px",
    }

    let off_style = {
        textAlign: "left",
        height:"80px",
        width:"200px",
    }

    const handleClick = (event) => {
        const new_val = !toggleOn;
        handleUpdate(new_val)

        if (new_val){
            setMyStyle(on_style)
            setMyText("On")
        } else {
            setMyStyle(off_style)
            setMyText("Off")
        }


    }


    const [myStyle, setMyStyle] = useState(toggleOn? on_style:off_style);
    const [myText, setMyText] = useState(toggleOn? "On":"Off");


    return (
        <div>
        <button style={myStyle} onClick={handleClick}>
            <p width="100%">{myText}</p>
        </button>
        </div>

    )
}
