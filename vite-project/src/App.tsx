import { useState } from 'react'
import reactLogo from './assets/react.svg'
import viteLogo from '/vite.svg'
import './App.css'
import SlideToggle from './components/slide_toggle/SlideToggle'

function App() {

  const [value, setValue] =useState(false);

    function handleUpdate(data){
      setValue(data);
      console.log(data)
    }

    const [slideValue, setSlideValue] =useState(false);

    function handleSlideToggleUpdate(data){
      setSlideValue(data);
      console.log(data)
    }

  return (
    <>
      <SlideToggle  isRound={false} handleUpdate={handleSlideToggleUpdate} toggleOn={slideValue}/>
    </>
  )
}

export default App
